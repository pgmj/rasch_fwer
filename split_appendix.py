"""Split the rendered manuscript docx into the main text and the appendix.

    python3 split_appendix.py index.docx

writes index-main.docx and index-appendix.docx beside the input. Render the
full docx first (with --profile blind for the anonymous version); the split is
done on the rendered file so that every cross-reference keeps the number it
was given in the full document ("Table S3" in the body, "Figure 1" and
"Equation 1" in the appendix).

The split point is the page break before the first level-1 "Appendix"
heading that apaquarto inserts ("Appendix" for one appendix, "Appendix A",
"Appendix B" and so on for several, or "Online Supplement" once
supplement.lua has relabelled it). Cross-reference links whose target ends up in the other
file are kept as plain text, bookmarks left without a partner are dropped,
and images the appendix does not use are removed from its package.

The reference list is split the same way. The main text keeps the entries it
cites, and the appendix gets its own list, on a new page at its end, with the
entries it cites. A work cited in both appears in both lists.
"""

import re
import sys
import zipfile
from pathlib import Path

PAGE_BREAK = '<w:p><w:r><w:br w:type="page"/></w:r></w:p>'
APPENDIX = re.compile(r'<w:p><w:pPr><w:pStyle w:val="Heading1" /></w:pPr><w:r>'
                      r'<w:t xml:space="preserve">(?:Appendix(?: [A-Z])?|Online Supplement)</w:t></w:r></w:p>')


# A reference entry is one Bibliography paragraph wrapped in its own bookmark,
# which the citations link to. Long citation keys get a hashed bookmark name,
# so entries are matched by bookmark, not by the "ref-" prefix.
ENTRY = re.compile(r'<w:bookmarkStart w:id="(\d+)" w:name="([^"]+)" />\s*'
                   r'<w:p><w:pPr><w:pStyle w:val="Bibliography" />'
                   r'(?:(?!</w:p>).)*</w:p>\s*<w:bookmarkEnd w:id="\1" />', re.S)
REFS_HEADING = ('<w:p><w:pPr><w:pStyle w:val="Heading1" /></w:pPr><w:r>'
                '<w:t xml:space="preserve">References</w:t></w:r></w:p>')


def cited(body):
    return set(re.findall(r'w:anchor="([^"]+)"', body))


def split_references(main_body, appx_body):
    """Keep main-text entries in the main list and give the appendix its own."""
    entries = list(ENTRY.finditer(main_body))
    if not entries:
        return main_body, appx_body
    first, last = entries[0].start(), entries[-1].end()
    in_main = cited(main_body[:first] + main_body[last:])
    in_appx = cited(appx_body)
    uncited = [e.group(2) for e in entries
               if e.group(2) not in in_main | in_appx]
    if uncited:
        sys.exit(f"Reference entries cited in neither part: {uncited}")
    main_list = "\n    ".join(e.group(0) for e in entries
                              if e.group(2) in in_main)
    appx_list = "\n    ".join(e.group(0) for e in entries
                              if e.group(2) in in_appx)
    main_body = main_body[:first] + main_list + main_body[last:]
    if appx_list:
        appx_body += ("\n    " + PAGE_BREAK + "\n    " + REFS_HEADING
                      + "\n    " + appx_list)
    return main_body, appx_body


def tidy(body):
    """Unlink dangling cross-references and drop unpaired bookmarks."""
    starts = set(re.findall(r'<w:bookmarkStart w:id="(\d+)"', body))
    ends = set(re.findall(r'<w:bookmarkEnd w:id="(\d+)"', body))
    body = re.sub(r'\s*<w:bookmarkStart w:id="(\d+)"[^>]*/>',
                  lambda m: m.group(0) if m.group(1) in ends else "", body)
    body = re.sub(r'\s*<w:bookmarkEnd w:id="(\d+)"\s*/>',
                  lambda m: m.group(0) if m.group(1) in starts else "", body)
    names = set(re.findall(r'<w:bookmarkStart [^>]*w:name="([^"]+)"', body))
    return re.sub(r'<w:hyperlink w:anchor="([^"]+)"[^>]*>(.*?)</w:hyperlink>',
                  lambda m: m.group(0) if m.group(1) in names else m.group(2),
                  body, flags=re.S)


def write(src, dst, document, drop_unused_media):
    rels_name = "word/_rels/document.xml.rels"
    with zipfile.ZipFile(src) as zin:
        rels = zin.read(rels_name).decode("utf8")
        dropped = set()
        if drop_unused_media:
            used = set(re.findall(r'r:(?:embed|id)="([^"]+)"', document))
            for m in re.finditer(r'<Relationship [^>]*/>', rels):
                rel = m.group(0)
                rid = re.search(r'Id="([^"]+)"', rel).group(1)
                target = re.search(r'Target="([^"]+)"', rel).group(1)
                if target.startswith("media/") and rid not in used:
                    dropped.add("word/" + target)
                    rels = rels.replace(rel, "")
        with zipfile.ZipFile(dst, "w", zipfile.ZIP_DEFLATED) as zout:
            for item in zin.infolist():
                if item.filename in dropped:
                    continue
                if item.filename == "word/document.xml":
                    data = document.encode("utf8")
                elif item.filename == rels_name:
                    data = rels.encode("utf8")
                else:
                    data = zin.read(item.filename)
                zout.writestr(item, data)


def main(path):
    src = Path(path)
    with zipfile.ZipFile(src) as z:
        xml = z.read("word/document.xml").decode("utf8")

    head, rest = xml.split("<w:body>", 1)
    body, tail = rest.rsplit("<w:sectPr", 1)
    tail = "<w:sectPr" + tail

    m = APPENDIX.search(body)
    if m is None:
        sys.exit("No level-1 'Appendix' or 'Online Supplement' heading found.")
    i = m.start()
    j = body.rfind(PAGE_BREAK, 0, i)
    if j < 0:
        sys.exit("No page break found before the Appendix heading.")

    # Bookmark ends between the page break and the heading close sections of
    # the main text, so they stay with it.
    main_body = body[:j] + body[j + len(PAGE_BREAK):i]
    appx_body = body[i:]
    main_body, appx_body = split_references(main_body, appx_body)

    out = {}
    for part, part_body in (("main", main_body), ("appendix", appx_body)):
        document = head + "<w:body>" + tidy(part_body) + "\n    " + tail
        dst = src.with_name(f"{src.stem}-{part}.docx")
        write(src, dst, document, drop_unused_media=(part == "appendix"))
        out[part] = dst
    for part, dst in out.items():
        print(f"{part:9s} {dst}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "index.docx")
