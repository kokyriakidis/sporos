"""Checks every DOI in refs.py against its registry: the title registered for it must match ours, and the year
must be the year of print (or of first posting for a preprint). Crossref holds journal and most proceedings DOIs;
DataCite holds arXiv's and LIPIcs's, and is asked when Crossref does not know a DOI.

    python3 tools/paper/check_refs.py
"""
import json, os, re, sys, urllib.request
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from refs import REFS

# Registered titles that differ from the printed one, by key: the registry's text, checked by hand.
REGISTERED = {
    "Schleimer2003": "Winnowing",  # Crossref keeps only the part before the colon
    "Hayashi2001": "Complete Genome Sequence of Enterohemorrhagic Eschelichia coli O157:H7 and Genomic Comparison "
                   "with a Laboratory Strain K-12",  # the journal's misprint, kept by Crossref
}
AGENT = {"User-Agent": "sporos-paper-check (mailto:kkyriaki@ucsc.edu)"}


def norm(s):
    """Letters and digits only, inline markup (<i>, <scp>) removed."""
    return re.sub(r"[^a-z0-9]", "", re.sub(r"<[^>]+>", "", s).lower())


def fetch(url):
    with urllib.request.urlopen(urllib.request.Request(url, headers=AGENT), timeout=30) as f:
        return json.load(f)


def registered(doi):
    """(registry, title, year, volume, pages) as Crossref or, failing it, DataCite holds them."""
    try:
        m = fetch(f"https://api.crossref.org/works/{doi}")["message"]
        year = (m.get("published-print") or m.get("issued") or {}).get("date-parts", [[None]])[0][0]
        return "crossref", (m.get("title") or [""])[0], year, m.get("volume"), m.get("page") or m.get("article-number") or ""
    except urllib.error.HTTPError as e:
        if e.code != 404:
            raise
    a = fetch(f"https://api.datacite.org/dois/{doi}")["data"]["attributes"]
    return "datacite", (a.get("titles") or [{}])[0].get("title", ""), a.get("publicationYear"), None, ""


bad = 0
for key, r in REFS.items():
    doi = r.get("doi")
    if not doi:
        continue
    try:
        where, title, year, volume, pages = registered(doi)
    except Exception as e:
        print(f"{key}: {doi}: {e}"); bad += 1; continue
    ok = norm(title) == norm(REGISTERED.get(key, r["title"]))
    print(f"{'ok ' if ok else 'BAD'} {key} ({where}): year {year} (ours {r['year']}), vol {volume} (ours {r.get('volume')}), "
          f"pages {pages} (ours {r.get('pages')})" + ("" if ok else f"\n    registered: {title}\n    ours:       {r['title']}"))
    bad += not ok
sys.exit(1 if bad else 0)
