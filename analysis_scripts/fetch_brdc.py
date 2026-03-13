import argparse
import datetime as dt
import gzip
import re
import shutil
import sys
from datetime import timezone
from pathlib import Path
from urllib.request import Request, urlopen

import certifi
import ssl

USER_AGENT = "rnx2rtkp-fetch/1.0"
SSL_CTX = ssl.create_default_context(cafile=certifi.where())


def http_get_bytes(url: str) -> bytes:
    req = Request(url, headers={"User-Agent": USER_AGENT})
    with urlopen(req, timeout=60, context=SSL_CTX) as resp:
        return resp.read()


def download_to(url: str, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    req = Request(url, headers={"User-Agent": USER_AGENT})
    with urlopen(req, timeout=60, context=SSL_CTX) as resp, open(out_path, "wb") as f:
        shutil.copyfileobj(resp, f)


def utc_date_today() -> dt.date:
    return dt.datetime.now(timezone.utc).date()


def doy(date: dt.date) -> int:
    return int(date.strftime("%j"))


def gunzip(gz_path: Path, out_path: Path) -> None:
    with gzip.open(gz_path, "rb") as fin, open(out_path, "wb") as fout:
        shutil.copyfileobj(fin, fout)


def pick_brdc_filename(index_html: str, year: int, day_of_year: int) -> str:
    tag = f"{year:04d}{day_of_year:03d}0000_01D_MN.rnx.gz"
    pattern = re.compile(rf'href="(BRDC00[^"]+_{tag})"', re.IGNORECASE)

    matches = pattern.findall(index_html)
    if not matches:
        return ""

    # Prefer "_R_" (rapid?) when present, otherwise take first match.
    for m in matches:
        if "_R_" in m:
            return m
    return matches[0]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--year", type=int, required=True)
    parser.add_argument("--month", type=int, required=True)
    parser.add_argument("--day", type=int, required=True)
    args = parser.parse_args()

    date = dt.date(args.year, args.month, args.day)
    year = date.year
    day_of_year = doy(date)

    base_dir = f"https://igs.bkg.bund.de/root_ftp/IGS/BRDC/{year:04d}/{day_of_year:03d}/"
    print("Index:", base_dir)

    index_html = http_get_bytes(base_dir).decode("utf-8", errors="replace")
    filename = pick_brdc_filename(index_html, year, day_of_year)
    if not filename:
        raise RuntimeError(
            f"Could not find BRDC00* file in index for {year}/{day_of_year:03d}"
        )

    file_url = base_dir + filename

    outdir = Path(args.outdir).resolve()
    gz_path = outdir / "brdc.rnx.gz"
    rnx_path = outdir / "brdc.rnx"

    print("Downloading:", file_url)
    download_to(file_url, gz_path)

    print("Decompressing:", gz_path, "->", rnx_path)
    gunzip(gz_path, rnx_path)

    print("Done:", rnx_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())