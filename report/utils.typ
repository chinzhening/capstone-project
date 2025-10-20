#let fonts = (
  serif: "Libertinus Serif",
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
  info: none,
) = {
  show: set par(spacing: 1.5em)
  v(1fr)
  align(center,{
    block(
      above: 2.5em,
      below: 2.5em,
      text(size: 16pt, font: fonts.sans, title)
    )

    // names
    block(text(size: 12pt, student))
    block(text(size: 12pt, supervisor))
    block(text(size: 12pt, examiner))
    
    // institution information
    block(text(size: 12pt, department))
    block(text(size: 12pt, institution))
    block(text(size: 12pt, info))
  })
  v(1fr)
}

