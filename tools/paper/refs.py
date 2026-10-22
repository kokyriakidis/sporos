"""The paper's references, written once: render.py makes paper/references.bib and the HTML reference list from
them. Every DOI is checked against Crossref by check_refs.py; entries without one are software or data.

Each entry: key -> dict(type, authors [family, given], title, venue fields, year, doi, cite (text citation)).
"""

def A(*names):
    """Authors as (family, given) pairs from 'Family, Given' strings; 'others' ends a truncated list."""
    return [tuple(n.split(", ", 1)) if ", " in n else (n, "") for n in names]


REFS = {
# Reference entries for the sporos paper. Names use Unicode; only "_" is LaTeX-escaped (Zheng2020 number).
# Verified 2026-10-21 against Crossref (journal/proceedings DOIs), DataCite (LIPIcs and arXiv DOIs)
# or Europe PMC (author lists where Crossref is truncated). Article numbers are given in `pages`.
# Entries without a DOI are software/tech reports; see the comments.

# ---------------------------------------------------------------- minimizers / sampling schemes
 "Roberts2004": dict(type="article", authors=A("Roberts, Michael", "Hayes, Wayne", "Hunt, Brian R.", "Mount, Stephen M.", "others"), title="Reducing storage requirements for biological sequence comparison", journal="Bioinformatics", year="2004", volume="20", number="18", pages="3363--3369", doi="10.1093/bioinformatics/bth408"),
 "Schleimer2003": dict(type="article", authors=A("Schleimer, Saul", "Wilkerson, Daniel S.", "Aiken, Alex"), title="Winnowing: local algorithms for document fingerprinting", journal="Proceedings of the 2003 ACM SIGMOD International Conference on Management of Data", year="2003", pages="76--85", doi="10.1145/872757.872770"),  # Crossref main title is just "Winnowing"; the subtitle is standard
 "Marcais2017": dict(type="article", authors=A("Marçais, Guillaume", "Pellow, David", "Bork, Daniel", "Orenstein, Yaron", "others"), title="Improving the performance of minimizers and winnowing schemes", journal="Bioinformatics", year="2017", volume="33", number="14", pages="i110--i117", doi="10.1093/bioinformatics/btx235"),
 "Zheng2020": dict(type="article", authors=A("Zheng, Hongyu", "Kingsford, Carl", "Marçais, Guillaume"), title="Improved design and analysis of practical minimizers", journal="Bioinformatics", year="2020", volume="36", number="Supplement\\_1", pages="i119--i127", doi="10.1093/bioinformatics/btaa472"),
 "GrootKoerkamp2024mod": dict(type="article", authors=A("Groot Koerkamp, Ragnar", "Pibiri, Giulio Ermanno"), title="The mod-minimizer: a simple and efficient sampling algorithm for long k-mers", journal="24th International Workshop on Algorithms in Bioinformatics (WABI 2024), LIPIcs", year="2024", volume="312", pages="11:1--11:23", doi="10.4230/LIPIcs.WABI.2024.11"),  # DataCite
 "GrootKoerkamp2025oc": dict(type="article", authors=A("Groot Koerkamp, Ragnar", "Liu, Daniel", "Pibiri, Giulio Ermanno"), title="The open-closed mod-minimizer algorithm", journal="Algorithms for Molecular Biology", year="2025", volume="20", number="1", pages="4", doi="10.1186/s13015-025-00270-0"),
 "GrootKoerkamp2025simd": dict(type="article", authors=A("Groot Koerkamp, Ragnar", "Martayan, Igor"), title="SimdMinimizers: computing random minimizers, fast", journal="23rd International Symposium on Experimental Algorithms (SEA 2025), LIPIcs", year="2025", volume="338", pages="20:1--20:19", doi="10.4230/LIPIcs.SEA.2025.20"),  # DataCite; full version bioRxiv 10.1101/2025.01.27.634998; no arXiv version found
 "GrootKoerkamp2025simdFull": dict(type="misc", authors=A("Groot Koerkamp, Ragnar", "Martayan, Igor"), title="SimdMinimizers: computing random minimizers, fast", journal="bioRxiv", year="2025", doi="10.1101/2025.01.27.634998"),  # optional: full version

# ---------------------------------------------------------------- syncmers / hashing
 "Edgar2021": dict(type="article", authors=A("Edgar, Robert"), title="Syncmers are more sensitive than minimizers for selecting conserved k-mers in biological sequences", journal="PeerJ", year="2021", volume="9", pages="e10805", doi="10.7717/peerj.10805"),
 "Shaw2022": dict(type="article", authors=A("Shaw, Jim", "Yu, Yun William"), title="Theory of local k-mer selection with applications to long-read alignment", journal="Bioinformatics", year="2022", volume="38", number="20", pages="4659--4669", doi="10.1093/bioinformatics/btab790"),  # online 2021, issue 2022
 "Mohamadi2016": dict(type="article", authors=A("Mohamadi, Hamid", "Chu, Justin", "Vandervalk, Benjamin P.", "Birol, Inanc"), title="ntHash: recursive nucleotide hashing", journal="Bioinformatics", year="2016", volume="32", number="22", pages="3492--3494", doi="10.1093/bioinformatics/btw397"),
 "Kazemi2022": dict(type="article", authors=A("Kazemi, Parham", "Wong, Johnathan", "Nikolić, Vladimir", "Mohamadi, Hamid", "others"), title="ntHash2: recursive spaced seed hashing for nucleotide sequences", journal="Bioinformatics", year="2022", volume="38", number="20", pages="4812--4813", doi="10.1093/bioinformatics/btac564"),

# ---------------------------------------------------------------- mappers using hash seed indexes
 "Li2018": dict(type="article", authors=A("Li, Heng"), title="Minimap2: pairwise alignment for nucleotide sequences", journal="Bioinformatics", year="2018", volume="34", number="18", pages="3094--3100", doi="10.1093/bioinformatics/bty191"),
 "Sahlin2021": dict(type="article", authors=A("Sahlin, Kristoffer"), title="Effective sequence similarity detection with strobemers", journal="Genome Research", year="2021", volume="31", number="11", pages="2080--2094", doi="10.1101/gr.275648.121"),
 "Sahlin2022": dict(type="article", authors=A("Sahlin, Kristoffer"), title="Strobealign: flexible seed size enables ultra-fast and accurate read alignment", journal="Genome Biology", year="2022", volume="23", number="1", pages="260", doi="10.1186/s13059-022-02831-7"),
 "Tolstoganov2026": dict(type="article", authors=A("Tolstoganov, Ivan", "Martin, Marcel", "Buchin, Nicolas", "Sahlin, Kristoffer"), title="Multi-context seeds enable fast and high-accuracy read mapping", journal="Genome Biology", year="2026", volume="27", number="1", pages="118", doi="10.1186/s13059-026-04017-x"),  # current strobealign citation
 "Karami2024": dict(type="article", authors=A("Karami, Moein", "Soltani Mohammadi, Aryan", "Martin, Marcel", "Ekim, Barış", "others"), title="Designing efficient randstrobes for sequence similarity analyses", journal="Bioinformatics", year="2024", volume="40", number="4", pages="btae187", doi="10.1093/bioinformatics/btae187"),  # optional
 "Rautiainen2020": dict(type="article", authors=A("Rautiainen, Mikko", "Marschall, Tobias"), title="GraphAligner: rapid and versatile sequence-to-graph alignment", journal="Genome Biology", year="2020", volume="21", number="1", pages="253", doi="10.1186/s13059-020-02157-2"),
 "Alser2021": dict(type="article", authors=A("Alser, Mohammed", "Rotman, Jeremy", "Deshpande, Dhrithi", "Taraszka, Kodi", "others"), title="Technology dictates algorithms: recent developments in read alignment", journal="Genome Biology", year="2021", volume="22", number="1", pages="249", doi="10.1186/s13059-021-02443-7"),  # optional survey

# ---------------------------------------------------------------- MEM / MUM finding between genomes
 "Delcher1999": dict(type="article", authors=A("Delcher, Arthur L.", "Kasif, Simon", "Fleischmann, Robert D.", "Peterson, Jeremy", "others"), title="Alignment of whole genomes", journal="Nucleic Acids Research", year="1999", volume="27", number="11", pages="2369--2376", doi="10.1093/nar/27.11.2369"),
 "Delcher2002": dict(type="article", authors=A("Delcher, Arthur L.", "Phillippy, Adam", "Carlton, Jane", "Salzberg, Steven L."), title="Fast algorithms for large-scale genome alignment and comparison", journal="Nucleic Acids Research", year="2002", volume="30", number="11", pages="2478--2483", doi="10.1093/nar/30.11.2478"),  # MUMmer2 (optional); Crossref lists only the first author, co-authors from the standard citation
 "Kurtz2004": dict(type="article", authors=A("Kurtz, Stefan", "Phillippy, Adam", "Delcher, Arthur L.", "Smoot, Michael", "others"), title="Versatile and open software for comparing large genomes", journal="Genome Biology", year="2004", volume="5", number="2", pages="R12", doi="10.1186/gb-2004-5-2-r12"),
 "Marcais2018": dict(type="article", authors=A("Marçais, Guillaume", "Delcher, Arthur L.", "Phillippy, Adam M.", "Coston, Rachel", "others"), title="MUMmer4: a fast and versatile genome alignment system", journal="PLOS Computational Biology", year="2018", volume="14", number="1", pages="e1005944", doi="10.1371/journal.pcbi.1005944"),
 "Abouelhoda2004": dict(type="article", authors=A("Abouelhoda, Mohamed Ibrahim", "Kurtz, Stefan", "Ohlebusch, Enno"), title="Replacing suffix trees with enhanced suffix arrays", journal="Journal of Discrete Algorithms", year="2004", volume="2", number="1", pages="53--86", doi="10.1016/S1570-8667(03)00065-0"),  # optional
 "Khan2009": dict(type="article", authors=A("Khan, Zia", "Bloom, Joshua S.", "Kruglyak, Leonid", "Singh, Mona"), title="A practical algorithm for finding maximal exact matches in large sequence datasets using sparse suffix arrays", journal="Bioinformatics", year="2009", volume="25", number="13", pages="1609--1616", doi="10.1093/bioinformatics/btp275"),
 "Vyverman2013": dict(type="article", authors=A("Vyverman, Michaël", "De Baets, Bernard", "Fack, Veerle", "Dawyndt, Peter"), title="essaMEM: finding maximal exact matches using enhanced sparse suffix arrays", journal="Bioinformatics", year="2013", volume="29", number="6", pages="802--804", doi="10.1093/bioinformatics/btt042"),
 "Khiste2015": dict(type="article", authors=A("Khiste, Nilesh", "Ilie, Lucian"), title="E-MEM: efficient computation of maximal exact matches for very large genomes", journal="Bioinformatics", year="2015", volume="31", number="4", pages="509--514", doi="10.1093/bioinformatics/btu687"),  # online 2014
 "Liu2019bfMEM": dict(type="article", authors=A("Liu, Yuansheng", "Zhang, Leo Yu", "Li, Jinyan"), title="Fast detection of maximal exact matches via fixed sampling of query K-mers and Bloom filtering of index K-mers", journal="Bioinformatics", year="2019", volume="35", number="22", pages="4560--4567", doi="10.1093/bioinformatics/btz273"),  # bfMEM (optional)
 "Grabowski2019": dict(type="article", authors=A("Grabowski, Szymon", "Bieniecki, Wojciech"), title="copMEM: finding maximal exact matches via sampling both genomes", journal="Bioinformatics", year="2019", volume="35", number="4", pages="677--678", doi="10.1093/bioinformatics/bty670"),  # online 2018
 "Grabowski2023": dict(type="article", authors=A("Grabowski, Szymon", "Bieniecki, Wojciech"), title="copMEM2: robust and scalable maximum exact match finding", journal="Bioinformatics", year="2023", volume="39", number="5", pages="btad313", doi="10.1093/bioinformatics/btad313"),
 "Shivakumar2025": dict(type="article", authors=A("Shivakumar, Vikram S.", "Langmead, Ben"), title="Mumemto: efficient maximal matching across pangenomes", journal="Genome Biology", year="2025", volume="26", number="1", pages="169", doi="10.1186/s13059-025-03644-0"),
 "Shivakumar2026": dict(type="article", authors=A("Shivakumar, Vikram S.", "Langmead, Ben"), title="Partitioned multi-MUM finding for scalable pangenomics with MumemtoM", journal="Genome Research", year="2026", volume="36", number="2", pages="397--404", doi="10.1101/gr.280940.125"),  # vol/issue/pages from Europe PMC (PMID 41203424); Crossref record lacks them; optional

# ---------------------------------------------------------------- BWT / FM-index / SMEM
 "BurrowsWheeler1994": dict(type="misc", authors=A("Burrows, Michael", "Wheeler, David J."), title="A block-sorting lossless data compression algorithm", journal="Technical Report 124, Digital Equipment Corporation Systems Research Center, Palo Alto, CA", year="1994"),  # no DOI; not in Crossref; standard citation
 "Ferragina2000": dict(type="article", authors=A("Ferragina, Paolo", "Manzini, Giovanni"), title="Opportunistic data structures with applications", journal="Proceedings of the 41st Annual Symposium on Foundations of Computer Science (FOCS)", year="2000", pages="390--398", doi="10.1109/SFCS.2000.892127"),
 "Navarro2007": dict(type="article", authors=A("Navarro, Gonzalo", "Mäkinen, Veli"), title="Compressed full-text indexes", journal="ACM Computing Surveys", year="2007", volume="39", number="1", pages="2", doi="10.1145/1216370.1216372"),  # optional
 "Lam2009": dict(type="article", authors=A("Lam, Tak-Wah", "Li, Ruiqiang", "Tam, Alan", "Wong, Simon", "others"), title="High throughput short read alignment via bi-directional BWT", journal="2009 IEEE International Conference on Bioinformatics and Biomedicine (BIBM)", year="2009", pages="31--36", doi="10.1109/BIBM.2009.42"),  # bidirectional BWT (optional); Crossref first author "Lam, T. W."
 "Li2009": dict(type="article", authors=A("Li, Heng", "Durbin, Richard"), title="Fast and accurate short read alignment with Burrows--Wheeler transform", journal="Bioinformatics", year="2009", volume="25", number="14", pages="1754--1760", doi="10.1093/bioinformatics/btp324"),
 "Li2012": dict(type="article", authors=A("Li, Heng"), title="Exploring single-sample SNP and INDEL calling with whole-genome de novo assembly", journal="Bioinformatics", year="2012", volume="28", number="14", pages="1838--1844", doi="10.1093/bioinformatics/bts280"),  # FMD-index / SMEM algorithm
 "Li2013": dict(type="misc", authors=A("Li, Heng"), title="Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM", journal="arXiv", year="2013", number="1303.3997", doi="10.48550/arXiv.1303.3997"),
 "Vasimuddin2019": dict(type="article", authors=A("Vasimuddin, Md", "Misra, Sanchit", "Li, Heng", "Aluru, Srinivas"), title="Efficient architecture-aware acceleration of BWA-MEM for multicore systems", journal="2019 IEEE International Parallel and Distributed Processing Symposium (IPDPS)", year="2019", pages="314--324", doi="10.1109/IPDPS.2019.00041"),  # arXiv:1907.12931
 "Misra2018": dict(type="article", authors=A("Misra, Sanchit", "Pan, Tony C.", "Mahadik, Kanak", "Powley, George", "others"), title="Performance extraction and suitability analysis of multi- and many-core architectures for next generation sequencing secondary analysis", journal="Proceedings of the 27th International Conference on Parallel Architectures and Compilation Techniques (PACT)", year="2018", pages="1--14", doi="10.1145/3243176.3243197"),  # batched + prefetched FM-index exact search (optional)
 "Li2024": dict(type="article", authors=A("Li, Heng"), title="BWT construction and search at the terabase scale", journal="Bioinformatics", year="2024", volume="40", number="12", pages="btae717", doi="10.1093/bioinformatics/btae717"),  # ropebwt3; arXiv:2409.00613
 "Li2026": dict(type="misc", authors=A("Li, Heng", "Homer, Nils"), title="Fast genomic read alignment with minibwa", journal="arXiv", year="2026", number="2606.15357", doi="10.48550/arXiv.2606.15357"),  # minibwa; v1 13 Jun 2026; code https://github.com/lh3/minibwa
 "Gagie2024": dict(type="article", authors=A("Gagie, Travis"), title="How to find long maximal exact matches and ignore short ones", journal="Developments in Language Theory (DLT 2024), Lecture Notes in Computer Science", year="2024", volume="14791", pages="131--140", doi="10.1007/978-3-031-66159-4_10"),  # LNCS volume from Springer series listing (not Crossref)
 "Gagie2026": dict(type="misc", authors=A("Gagie, Travis"), title="Parse indexing for discarding short pseudo-MEMs safely", journal="arXiv", year="2026", number="2605.17574", doi="10.48550/arXiv.2605.17574"),  # defines f-MEMs in its Introduction (v1 17 May 2026, v3 28 May 2026); not a paper titled after f-MEMs
 "Tatarnikov2022": dict(type="misc", authors=A("Tatarnikov, Igor", "Shahrabi Farahani, Ardavan", "Kashgouli, Sana", "Gagie, Travis"), title="MONI can find k-MEMs", journal="arXiv", year="2022", number="2202.05085", doi="10.48550/arXiv.2202.05085"),  # earlier name for f-MEMs; no peer-reviewed version found
 "Brown2025": dict(type="article", authors=A("Brown, Nathaniel K.", "Depuydt, Lore", "Zakeri, Mohsen", "Alhadi, Anas", "others"), title="KeBaB: k-mer based breaking for finding long MEMs", journal="String Processing and Information Retrieval (SPIRE 2025), Lecture Notes in Computer Science", year="2025", pages="10--17", doi="10.1007/978-3-032-05228-5_2"),  # UNVERIFIED: LNCS volume number not found; DOI/pages/authors verified (optional)
 "Rossi2022": dict(type="article", authors=A("Rossi, Massimiliano", "Oliva, Marco", "Langmead, Ben", "Gagie, Travis", "others"), title="MONI: a pangenomic index for finding maximal exact matches", journal="Journal of Computational Biology", year="2022", volume="29", number="2", pages="169--187", doi="10.1089/cmb.2021.0290"),
 "Gagie2020": dict(type="article", authors=A("Gagie, Travis", "Navarro, Gonzalo", "Prezza, Nicola"), title="Fully functional suffix trees and optimal text searching in BWT-runs bounded space", journal="Journal of the ACM", year="2020", volume="67", number="1", pages="2:1--2:54", doi="10.1145/3375890"),  # Crossref pages 1-54 (article 2)

# ---------------------------------------------------------------- suffix-array construction
 "Nong2009": dict(type="article", authors=A("Nong, Ge", "Zhang, Sen", "Chan, Wai Hong"), title="Linear suffix array construction by almost pure induced-sorting", journal="2009 Data Compression Conference (DCC)", year="2009", pages="193--202", doi="10.1109/DCC.2009.42"),
 "Nong2011": dict(type="article", authors=A("Nong, Ge", "Zhang, Sen", "Chan, Wai Hong"), title="Two efficient algorithms for linear time suffix array construction", journal="IEEE Transactions on Computers", year="2011", volume="60", number="10", pages="1471--1484", doi="10.1109/TC.2010.188"),
 "libsais": dict(type="misc", authors=A("Grebnov, Ilya"), title="libsais: fast linear-time construction of suffix array, LCP array and Burrows--Wheeler transform", year="2026", url="https://github.com/IlyaGrebnov/libsais"),  # software, no paper; latest release v2.12.1 (2026-10-21) -- set year/version to the one used

# ---------------------------------------------------------------- seeding-share evidence
 "Subramaniyan2021": dict(type="article", authors=A("Subramaniyan, Arun", "Wadden, Jack", "Goliya, Kush", "Ozog, Nathan", "others"), title="Accelerated seeding for genome sequence alignment with enumerated radix trees", journal="2021 ACM/IEEE 48th Annual International Symposium on Computer Architecture (ISCA)", year="2021", pages="388--401", doi="10.1109/ISCA52012.2021.00038"),  # "~40% of BWA-MEM2 time" (abstract, full text not accessed)
 "Kalikar2022": dict(type="article", authors=A("Kalikar, Saurabh", "Jain, Chirag", "Vasimuddin, Md", "Misra, Sanchit"), title="Accelerating minimap2 for long-read sequencing applications on modern CPUs", journal="Nature Computational Science", year="2022", volume="2", number="2", pages="78--83", doi="10.1038/s43588-022-00201-8"),  # mm2-fast; seeding 3-13% (bioRxiv text)
 "Eudine2023": dict(type="misc", authors=A("Eudine, Julien", "Alser, Mohammed", "Singh, Gagandeep", "Alkan, Can", "others"), title="GateSeeder: near-memory CPU-FPGA acceleration of short and long read mapping", journal="arXiv", year="2023", number="2309.17063", doi="10.48550/arXiv.2309.17063"),  # UNVERIFIED: peer-reviewed version (if any) not found; arXiv metadata verified

# ---------------------------------------------------------------- language / compiler
 "Mojo": dict(type="misc", authors=A("Modular Inc."), title="Mojo: a programming language for heterogeneous compute", year="2026", url="https://www.modular.com/mojo"),  # software, no paper; docs https://docs.modular.com/mojo/ ; source https://github.com/modular/modular -- set year/version to the one used; title wording is ours
 "Lattner2021": dict(type="article", authors=A("Lattner, Chris", "Amini, Mehdi", "Bondhugula, Uday", "Cohen, Albert", "others"), title="MLIR: scaling compiler infrastructure for domain specific computation", journal="2021 IEEE/ACM International Symposium on Code Generation and Optimization (CGO)", year="2021", pages="2--14", doi="10.1109/CGO51591.2021.9370308"),

# ---------------------------------------------------------------- data
 "Nurk2022": dict(type="article", authors=A("Nurk, Sergey", "Koren, Sergey", "Rhie, Arang", "Rautiainen, Mikko", "others"), title="The complete sequence of a human genome", journal="Science", year="2022", volume="376", number="6588", pages="44--53", doi="10.1126/science.abj6987"),
 "Schneider2017": dict(type="article", authors=A("Schneider, Valerie A.", "Graves-Lindsay, Tina", "Howe, Kerstin", "Bouk, Nathan", "others"), title="Evaluation of GRCh38 and de novo haploid genome assemblies demonstrates the enduring quality of the reference assembly", journal="Genome Research", year="2017", volume="27", number="5", pages="849--864", doi="10.1101/gr.213611.116"),
 "Yoo2025": dict(type="article", authors=A("Yoo, DongAhn", "Rhie, Arang", "Hebbar, Prajna", "Antonacci, Francesca", "others"), title="Complete sequencing of ape genomes", journal="Nature", year="2025", volume="641", number="8062", pages="401--418", doi="10.1038/s41586-025-08816-3"),
 "Blattner1997": dict(type="article", authors=A("Blattner, Frederick R.", "Plunkett, Guy", "Bloch, Craig A.", "Perna, Nicole T.", "others"), title="The complete genome sequence of Escherichia coli K-12", journal="Science", year="1997", volume="277", number="5331", pages="1453--1462", doi="10.1126/science.277.5331.1453"),
 "Hayashi2001": dict(type="article", authors=A("Hayashi, T.", "Makino, K.", "Ohnishi, M.", "Kurokawa, K.", "others"), title="Complete genome sequence of enterohemorrhagic Escherichia coli O157:H7 and genomic comparison with a laboratory strain K-12", journal="DNA Research", year="2001", volume="8", number="1", pages="11--22", doi="10.1093/dnares/8.1.11"),  # authors (initials only) from Europe PMC, PMID 11258796; Crossref lists only the first author; the journal record spells it "Eschelichia" (typo), title here is corrected
}


def cite(key):
    """(names, year) as an author-year citation shows them: one or two families, else 'First et al.'."""
    r = REFS[key]
    fam = [f for f, _ in r["authors"] if f != "others"]
    many = len(r["authors"]) > 2 or any(f == "others" for f, _ in r["authors"])
    names = fam[0] if len(fam) == 1 and not many else (f"{fam[0]} and {fam[1]}" if not many else f"{fam[0]} et al.")
    return names, r["year"]


def bibtex():
    out = ["% Generated by tools/paper/render.py from refs.py; edit the source there, not this file.",
           "% Every DOI is checked against Crossref by tools/paper/check_refs.py.", ""]
    for key, r in REFS.items():
        def esc(s):
            return s.replace("&", "\\&")
        authors = " and ".join("others" if f == "others" else (f"{{{f}}}" if not g else f"{f}, {g}") for f, g in r["authors"])
        fields = [("author", authors), ("title", "{" + esc(r["title"]) + "}")]
        if r["type"] == "article":
            fields.append(("journal", esc(r["journal"])))
        elif r.get("journal"):
            fields.append(("howpublished", esc(r["journal"])))
        for f in ("year", "volume", "number", "pages", "doi", "url"):
            if r.get(f):
                fields.append((f, r[f]))
        body = ",\n".join(f"  {k:<12} = {{{v}}}" for k, v in fields)
        out.append(f"@{r['type']}{{{key},\n{body}\n}}\n")
    return "\n".join(out)


def html_entry(key):
    """Authors (year) Title. Venue, volume, pages. doi."""
    r = REFS[key]
    names = [("et al." if f == "others" else (f"{f}," + "".join(p[0] + "." for p in g.replace("-", " ").split()) if g else f))
             for f, g in r["authors"]]
    if len(names) > 1 and names[-1] == "et al.":
        authors = ", ".join(names[:-1]) + " et al."
    else:
        authors = ", ".join(names[:-1]) + (" and " if len(names) > 1 else "") + names[-1]
    venue = r.get("journal", "")
    for f, sep in (("volume", ", "), ("pages", ", ")):
        if r.get(f):
            venue += sep + r[f].replace("--", "–")
    return authors, r["title"] + ".", venue, r.get("doi"), r.get("url")
