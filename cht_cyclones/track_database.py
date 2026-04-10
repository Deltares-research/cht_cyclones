"""
Cyclone track database: top-level container that holds one or more CycloneTrackDataset
instances and optionally synchronises them against a remote S3 bucket.
"""

import os

import boto3
import toml
from botocore import UNSIGNED
from botocore.client import Config

from cht_cyclones.track_dataset import CycloneTrackDataset


class CycloneTrackDatabase:
    """
    Top-level container for one or more tropical-cyclone track datasets.

    Reads dataset metadata from a local TOML catalogue file and optionally
    synchronises new datasets from a remote S3 bucket.

    Parameters
    ----------
    path : str
        Local directory where dataset folders and the catalogue file reside.
    s3_bucket : str, optional
        Name of the S3 bucket used for remote synchronisation.
    s3_key : str, optional
        Key prefix (folder path) inside the S3 bucket.
    s3_region : str, optional
        AWS region of the S3 bucket.
    check_online : bool, optional
        When ``True``, query the remote S3 bucket for new datasets on
        initialisation.
    """

    def __init__(
        self,
        path: str,
        s3_bucket: str | None = None,
        s3_key: str | None = None,
        s3_region: str | None = None,
        check_online: bool = False,
    ) -> None:
        self.path = path
        self.dataset = []
        self.s3_client = None
        self.s3_bucket = s3_bucket
        self.s3_key = s3_key
        self.s3_region = s3_region
        self.read()
        if check_online:
            self.check_online_database()

    def read(self) -> None:
        """
        Read metadata for all datasets listed in the local catalogue file.

        The catalogue file is expected at ``<path>/cyclone_track_datasets.tml``.
        Missing or unreadable dataset metadata files are skipped with a warning.
        """
        # Check if the path exists. If not, create it.
        if not os.path.exists(self.path):
            os.makedirs(self.path)

        # Read in database
        tml_file = os.path.join(self.path, "cyclone_track_datasets.tml")
        if not os.path.exists(tml_file):
            print(f"Warning! Cyclone tracks database file not found: {tml_file}")
            return

        datasets = toml.load(tml_file)

        for d in datasets["dataset"]:
            name = d["name"]

            if "path" in d:
                path = d["path"]
            else:
                path = os.path.join(self.path, name)

            # Read the meta data for this dataset
            fname = os.path.join(path, "metadata.tml")

            if os.path.exists(fname):
                metadata = toml.load(fname)
                dataset_format = metadata["format"]
            else:
                print(
                    f"Could not find metadata file for dataset {name} ! Skipping dataset."
                )
                continue

            if dataset_format.lower() == "ibtracs":
                dataset = CycloneTrackDataset(name, path)
            elif dataset_format.lower() == "hurdat2":
                pass

            self.dataset.append(dataset)

    def check_online_database(self) -> None:
        """
        Synchronise the local database with the remote S3 catalogue.

        Downloads the remote catalogue file, compares it with the local list of
        datasets, and fetches metadata for any new datasets found on S3.  The
        local catalogue TOML is updated if new datasets were added.
        """
        if self.s3_client is None:
            self.s3_client = boto3.client(
                "s3", config=Config(signature_version=UNSIGNED)
            )
        if self.s3_bucket is None:
            return
        # First download a copy of cyclone_track_datasets.tml and call it cyclone_track_datasets_s3.tml
        key = f"{self.s3_key}/cyclone_track_datasets.tml"
        filename = os.path.join(self.path, "cyclone_track_datasets_s3.tml")
        print("Updating cyclone track database ...")
        try:
            self.s3_client.download_file(
                Bucket=self.s3_bucket,  # assign bucket name
                Key=key,  # key is the file name
                Filename=filename,
            )  # storage file path
        except Exception:
            # Download failed
            print(
                f"Failed to download {key} from {self.s3_bucket}. Database will not be updated."
            )
            return

        # Read bathymetry_s3.tml
        short_name_list, long_name_list = self.dataset_names()
        datasets_s3 = toml.load(filename)
        track_datasets_added = False
        added_names = []
        # Loop through s3 datasets, and check whether they exist in the local database.
        # If so, check if the metadata also exists. If not, make local folder and download the metadata.
        # Additionally, check if available_tiles.nc in s3 and not in local database, download it.
        for d in datasets_s3["dataset"]:
            # Get list of existing datasets
            s3_name = d["name"]
            if s3_name not in short_name_list:
                # Dataset not in local database
                print(f"Adding track dataset {s3_name} to local database ...")
                # Create folder and download metadata
                path = os.path.join(self.path, s3_name)
                os.makedirs(path, exist_ok=True)
                key = f"{self.s3_key}/{s3_name}/metadata.tml"
                filename = os.path.join(path, "metadata.tml")
                # Download metadata
                try:
                    self.s3_client.download_file(
                        Bucket=self.s3_bucket,  # assign bucket name
                        Key=key,  # key is the file name
                        Filename=filename,
                    )  # storage file path
                except Exception as e:
                    print(e)
                    print(f"Failed to download {key}. Skipping tide model.")
                    continue
                # Necessary data has been downloaded
                track_datasets_added = True
                added_names.append(s3_name)
        # Write new local bathymetry.tml
        if track_datasets_added:
            d = {}
            d["dataset"] = []
            for name in short_name_list:
                d["dataset"].append({"name": name})
            for name in added_names:
                d["dataset"].append({"name": name})
            # Now write the new bathymetry.tml
            with open(
                os.path.join(self.path, "cyclone_track_datasets.tml"), "w"
            ) as tml:
                toml.dump(d, tml)
            # Read the database again
            self.dataset = []
            self.read()
        # else:
        #     print("No new tide models were added to the local database.")

    def get_dataset(self, name: str) -> "CycloneTrackDataset | None":
        """
        Retrieve a dataset by short name, downloading it first if necessary.

        Parameters
        ----------
        name : str
            Short name of the dataset as recorded in the catalogue.

        Returns
        -------
        CycloneTrackDataset or None
            The requested dataset, or ``None`` if the name was not found.
        """
        for dataset in self.dataset:
            if dataset.name == name:
                # Make sure the dataset is locally available
                dataset.download()
                dataset.read()
                return dataset
        return None

    def dataset_names(self) -> tuple[list[str], list[str]]:
        """
        Return the short and long names of all registered datasets.

        Returns
        -------
        short_name_list : list of str
            Short (identifier) names.
        long_name_list : list of str
            Human-readable long names.
        """
        short_name_list = []
        long_name_list = []
        for dataset in self.dataset:
            short_name_list.append(dataset.name)
            long_name_list.append(dataset.long_name)
        return short_name_list, long_name_list
