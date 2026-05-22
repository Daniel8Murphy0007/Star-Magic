-- Lua filter to convert \( and \) back to $ and $ for pdflatex compatibility
-- This filter returns Math elements as-is to let pandoc's default latex writer handle them correctly

function Math(elem)
  -- For InlineMath, return as inline TeX that will render correctly
  if elem.mathtype == 'InlineMath' then
    -- Return as a string that pdflatex can parse
    return pandoc.Str('') -- Let pandoc handle the math element normally
  end
  -- For DisplayMath, handle similarly
  return elem
end

-- Alternative approach: use a RawBlock filter to post-process the latex output
function RawBlock(elem)
  if elem.format == 'latex' then
    -- Replace \( with $ and \) with $
    local text = elem.text:gsub('\\(', '$'):gsub('\\)', '$')
    return pandoc.RawBlock('latex', text)
  end
  return elem
end

function RawInline(elem)
  if elem.format == 'latex' then
    -- Replace \( with $ and \) with $
    local text = elem.text:gsub('\\(', '$'):gsub('\\)', '$')
    return pandoc.RawInline('latex', text)
  end
  return elem
end
