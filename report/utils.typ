#let fonts = (
  serif: "CMU Serif",
  sans: "Noto Sans",
  mono: "Incosolata",
)

#let cover_page(
  title: none,
  student: none,
  supervisor: none,
  examiner: none,
  department: none,
  institution: none,
  module: none,
  info: none,
) = {
  show: set par(spacing: 1.5em)
  v(0.15fr)

  align(center,{
    // institution name
    block(text(size: 18pt, smallcaps[#institution]))

    v(0.15fr)

    // module code
    block(text(size: 16pt, smallcaps[#module]))

    v(0.1fr)

    line(length: 100%)
    block(
      above: 1.5em,
      below: 2.5em,
      text(size: 20pt, weight:"bold", title)
    )
    line(length: 100%)

    // names
    table(
      columns: (auto, 1fr, auto),
      inset: (x: 2em, y: 0.5em),
      stroke: none,
      align(left)[_Student_], [], align(right)[_Supervisor_],
      block(
        align(left,
        text(size: 12pt, smallcaps(student)))
      ),
      [],
      block(
        align(right,
        text(size: 12pt, weight: "thin", smallcaps(supervisor))
        )
      )
    )
    
    v(0.6fr)
    
    block(text(size: 12pt, examiner))
   
    block(text(size: 12pt, info))
  })
  v(0.3fr)
}