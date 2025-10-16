import os
import re
import hashlib
import requests
import feedparser
from bs4 import BeautifulSoup

RSS_URL = "https://www.metoc.navy.mil/jtwc/rss/jtwc.rss"
HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
        "AppleWebKit/537.36 (KHTML, like Gecko) "
        "Chrome/115.0.0.0 Safari/537.36"
    )
}
DOWNLOAD_DIR = "jmv3_downloads"
SEEN_FILE = "seen_jmv3.txt"

def load_seen():
    if os.path.exists(SEEN_FILE):
        with open(SEEN_FILE) as f:
            return set(line.strip() for line in f)
    return set()

def save_seen(seen):
    with open(SEEN_FILE, "w") as f:
        f.write("\n".join(sorted(seen)))

def md5(s):
    return hashlib.md5(s.encode("utf-8")).hexdigest()

def download_file(url, folder=DOWNLOAD_DIR):
    os.makedirs(folder, exist_ok=True)
    filename = os.path.join(folder, url.split("/")[-1])
    r = requests.get(url, headers=HEADERS)
    r.raise_for_status()
    with open(filename, "wb") as f:
        f.write(r.content)
    print(f"Downloaded: {filename}")

def main():
    seen = load_seen()
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
                if uid not in seen:
                    print(f"Found new JMV3.0: {href}")
                    try:
                        download_file(href)
                        seen.add(uid)
                        new_urls.append(href)
                    except Exception as e:
                        print(f"Error downloading {href}: {e}")

    if not new_urls:
        print("No new JMV3.0 files found.")
    save_seen(seen)

def download_jtwc_jmv30(path):
    """Download latest JTWc JMV 3.0 files from RSS feed.
    """
    DOWNLOAD_DIR = path
    main()

if __name__ == "__main__":
    main()
