--- Apply apaquarto's FigureNote paragraph style to table and figure notes.
---
--- Notes are written in the source as ordinary paragraphs beginning with an
--- emphasised "Note.", not as custom-style divs. A div would be the obvious
--- way to style them, but Quarto emits a div as a Typst `#block`, and a block
--- ends the paragraph run, so the paragraph after every note loses its
--- first-line indent in the PDF. Setting `all: true` on the Typst paragraph
--- indent fixes that but wrongly indents the text continuing a sentence after
--- a display equation.
---
--- Doing it here instead keeps one source form that is correct in Typst and
--- adds the Word styling only where it is wanted.

if not FORMAT:match("docx") then
  return {}
end

local function is_note(para)
  local first = para.content[1]
  return first ~= nil
    and first.t == "Emph"
    and pandoc.utils.stringify(first) == "Note."
end

function Para(el)
  if is_note(el) then
    return pandoc.Div({ el }, pandoc.Attr("", {}, { ["custom-style"] = "FigureNote" }))
  end
end
