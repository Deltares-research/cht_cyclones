import os
import shutil
import re
import hashlib
import requests
import feedparser
from bs4 import BeautifulSoup

from cht_cyclones import TropicalCycloneTrack, TropicalCyclone

RSS_URL = "https://www.metoc.navy.mil/jtwc/rss/jtwc.rss"
HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
        "AppleWebKit/537.36 (KHTML, like Gecko) "
        "Chrome/115.0.0.0 Safari/537.36"
    )
}
# DOWNLOAD_DIR = "jmv3_downloads"
# SEEN_FILE = "seen_jmv3.txt"

# def load_seen():
#     if os.path.exists(SEEN_FILE):
#         with open(SEEN_FILE) as f:
#             return set(line.strip() for line in f)
#     return set()

# def save_seen(seen):
#     with open(SEEN_FILE, "w") as f:
#         f.write("\n".join(sorted(seen)))

def md5(s):
    return hashlib.md5(s.encode("utf-8")).hexdigest()

def download_file(url, folder):
    os.makedirs(folder, exist_ok=True)
    filename = os.path.join(folder, url.split("/")[-1])
    r = requests.get(url, headers=HEADERS)
    r.raise_for_status()
    with open(filename, "wb") as f:
        f.write(r.content)
    print(f"Downloaded: {filename}")

def download_jmv30_files(path):

    # If path does not exist, create it
    if not os.path.exists(path):
        os.makedirs(path)

    r = requests.get(RSS_URL, headers=HEADERS)
    r.raise_for_status()
    
    feed = feedparser.parse(r.content)
    if feed.bozo:
        raise RuntimeError("RSS parsing failed")

    new_urls = []
    for entry in feed.entries:
        if "description" not in entry:
            continue
        
        # parse HTML inside <description>
        soup = BeautifulSoup(entry.description, "html.parser")
        for a in soup.find_all("a", href=True):
            href = a["href"]
            if href.lower().endswith(".tcw"):   # JMV 3.0 files
                uid = md5(href)
                print(f"Found new JMV3.0: {href}")
                try:
                    download_file(href, path)
                    new_urls.append(href)
                except Exception as e:
                    print(f"Error downloading {href}: {e}")

    if not new_urls:
        print("No new JMV3.0 files found.")

def organize(jtwc_path, download_path):

    # Loop through all downloaded files

    for filename in os.listdir(download_path):

        if filename.endswith(".tcw"):

            # Move or process the file as needed

            basin = filename[0:2]
            storm_num = filename[2:4]
            year = "20" + filename[4:6]

            if int(storm_num) >= 90:
                # This is not (yet) a named storm. Continue to next file.
                continue 

            track = TropicalCycloneTrack()
            config, name, advisory = track.read(os.path.join(download_path, filename), format="jmv30")
            advstr = f"{advisory:02d}" if advisory is not None else "xx"

            # Copy raw data file to jmv30 folder
            pth = os.path.join(jtwc_path, year, basin, storm_num, "jmv30")
            os.makedirs(pth, exist_ok=True)
            # Copy downloaded file to organized folder
            shutil.copy(os.path.join(download_path, filename), pth)            
            
            pth = os.path.join(jtwc_path, year, basin, storm_num)
            os.makedirs(pth, exist_ok=True)
            # And write *.cyc file
            track.write(os.path.join(pth, f"{basin}_{year}_{storm_num}_adv{advstr}.cyc"))

            # And now merge
            # Make a list of all existing cyc files in pth, include the complete path
            cyc_files = [os.path.join(pth, f) for f in os.listdir(pth) if f.endswith(".cyc")]
            tc = TropicalCyclone(track_file=cyc_files)
            tc.track.write(os.path.join(pth, f"{basin}_{year}_{storm_num}_merged.cyc"))
