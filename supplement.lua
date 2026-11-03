--- Present apaquarto's appendices as one online supplement in the docx.
---
--- The journal takes the appendix as a separate online supplement, so its
--- tables are numbered S1, S2 and so on rather than A1, A2. apaquarto gives
--- every level-1 section after the references its own appendix letter, which
--- it stores as a `prefix` attribute on each table and heads with an
--- "Appendix A", "Appendix B" label. Running after apaquarto's pre-ast
--- filters, this sets every table prefix to "S", relabels the first appendix
--- "Online Supplement", and drops the later labels. Cross-references and
--- captions are built from the prefix by apaquarto's later filters, so they
--- follow. This assumes only the first appendix holds tables, since the
--- table count restarts at every appendix.

if not FORMAT:match("docx") then
  return {}
end

local seen = false

local function relabel(h)
  if h.level == 1 and h.identifier == ""
      and pandoc.utils.stringify(h.content):match("^Appendix %u$") then
    if seen then
      return {}
    end
    seen = true
    h.content = pandoc.Inlines({ pandoc.Str("Online"), pandoc.Space(), pandoc.Str("Supplement") })
    return h
  end
end

local function reprefix(el)
  if el.identifier and el.identifier:match("^tbl%-") and el.attributes.prefix
      and el.attributes.prefix ~= "" then
    el.attributes.prefix = "S"
    return el
  end
end

return {
  {
    Header = relabel,
    FloatRefTarget = reprefix,
    Div = reprefix,
  }
}
