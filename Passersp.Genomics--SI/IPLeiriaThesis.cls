% Author Name: José António Portela Areia 
% Author Contact: jose.apareia@gmail.com
% Version: 2.0.0 - 2024/08/08
% Public Repository: https://github.com/joseareia/ipleiria-thesis

% Class information
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{IPLeiriaThesis}[2023/12/12 IPLeiria Thesis Class]

\hyphenpenalty=10000  % Prevent hyphenation
\exhyphenpenalty=10000 % Prevent hyphenation after explicit hyphens

% Process class options: document language
\DeclareOption{en}{\def\IPLeiria@language{english}}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{report}}
\ProcessOptions\relax

% Load master class
\LoadClass[a4paper,twoside,11pt]{report}

% Load the babel package with language settings
\RequirePackage[main=english]{babel}

\usepackage{polyglossia}
\setdefaultlanguage{english} % Primary document language
\setotherlanguage{portuguese} % Additional language

% Define Portuguese month names
\newcommand{\mespor}[1]{%
  \ifcase#1\or
    janeiro\or fevereiro\or março\or abril\or maio\or junho\or
    julho\or agosto\or setembro\or outubro\or novembro\or dezembro\fi}

% Define a command to get the selected language
\newcommand{\getLanguage}{\IPLeiria@language}
\ifx\IPLeiria@language\@undefined\renewcommand{\getLanguage}{english}\fi

% Select the language specified by the user
\AtBeginDocument{%
    \ifx\IPLeiria@language\@undefined
        \PackageWarning{IPLeiriaThesis}{Language not specified. Defaulting to English.}%
        \selectlanguage{portuguese}%
    \else
        \selectlanguage{\IPLeiria@language}%
    \fi
}

% Required packages
\usepackage{graphicx} % Required for including images
\usepackage{amsmath} % Mathematical expressions
\usepackage{amssymb} % Mathematical symbols
\usepackage[dvipsnames]{xcolor} % Using colors within the document
\usepackage{tocbibind} % Add list of figures and tables of contents
\usepackage[toc,page]{appendix} % Add Appendix
\usepackage{titling} % Titlepage and titleback customisation
\usepackage[calc]{datetime2} % Current date with both 
\usepackage{tcolorbox} % Highlighted box
\usepackage{subcaption} % For side-by-side figures
\usepackage{xparse} % "Davide" problem... solved!
\usepackage{ifthen} % Conditional statements
\usepackage{blindtext} % Generate dummy text
\usepackage{silence} % No warnings in certain packages
\usepackage{typearea} % KOMA options for page rotation
\usepackage{eso-pic} % For the background in both front and back matter
\usepackage{setspace} % For line height adjusts
\usepackage{varwidth} % Dependency for the custom boxes
\usepackage{pdfpages} % Allows the inclusion of whole PDFs
\usepackage{fontawesome5}
\usepackage{academicons}
\definecolor{orcidlogocol}{HTML}{A6CE39}

% Load more libaries from tcolorbox package
\tcbuselibrary{most}

% Silence some warnings
\WarningsOff[blindtext] % No support for portuguese language

% Colors
\colorlet{darkred}{red!50!black}
\definecolor{frontpagedark}{HTML}{1c1c1c}

% Margins
\usepackage[
	top=2cm, % Top margin
	bottom=2.5cm, % Bottom margin
	left=2.5cm, % Left margin
	right=2.5cm, % Right margin
	footskip=1.6cm, % Space from the bottom margin to footer
	headsep=0.75cm, % Space from the top margin to header
	% showframe % Uncomment to show frames around the margins for debugging
]{geometry}

% Fonts
% \usepackage[utf8]{inputenc} % Required for international characters
\usepackage[T1]{fontenc} % Font encoding for international characters
% \usepackage{newpxtext} % Alternative use of the PX fonts
\usepackage{newpxmath} % Alternative use of the PX fonts (Math)
% \usepackage{microtype} % Slightly tweak font spacing for aesthetics
\usepackage{fvextra} % Compatibility with csquotes (?)
\usepackage{csquotes} % Compatibility with babel (?)
\linespread{1.2} % Increase line spacing slightly

% Ensure compilation with XeLaTeX or LuaLaTeX for fontspec
\usepackage{fontspec}
\setmainfont{Cormorant Garamond}[LetterSpace=0.7]

% Defines a command to switch to Lato font
\newcommand{\latofont}{\fontspec{Cormorant Garamond}}

% Defines ORCID command ~
\newcommand{\orcid}[1]{\href{https://orcid.org/#1}{\textcolor[HTML]{A6CE39}{\aiOrcid}}}

% Document variables
\NewDocumentCommand{\firstauthor}{m}{\newcommand{\firstauthorname}{#1}}
\NewDocumentCommand{\firstauthorid}{m}{\newcommand{\firstauthornum}{#1}}

\NewDocumentCommand{\secondauthor}{m}{\newcommand{\secondauthorname}{#1}}
\NewDocumentCommand{\secondauthorid}{m}{\newcommand{\secondauthornum}{#1}}

\NewDocumentCommand{\thirdauthor}{m}{\newcommand{\thirdauthorname}{#1}}
\NewDocumentCommand{\thirdauthorid}{m}{\newcommand{\thirdauthornum}{#1}}

\NewDocumentCommand{\supervisor}{m}{\newcommand{\supname}{#1}}
\NewDocumentCommand{\supervisormail}{m}{\newcommand{\supmail}{#1}}
\NewDocumentCommand{\supervisortitle}{m}{\newcommand{\suptitle}{#1}}

\NewDocumentCommand{\cosupervisor}{m}{\newcommand{\cosupname}{#1}}
\NewDocumentCommand{\cosupervisormail}{m}{\newcommand{\cosupmail}{#1}}
\NewDocumentCommand{\cosupervisortitle}{m}{\newcommand{\cosuptitle}{#1}}

\NewDocumentCommand{\seccosupervisor}{m}{\newcommand{\seccosupname}{#1}}
\NewDocumentCommand{\seccosupervisormail}{m}{\newcommand{\seccosupmail}{#1}}
\NewDocumentCommand{\seccosupervisortitle}{m}{\newcommand{\seccosuptitle}{#1}}

\NewDocumentCommand{\subtitle}{m}{\newcommand{\subname}{#1}}
\NewDocumentCommand{\university}{m}{\newcommand{\univname}{#1}}
\NewDocumentCommand{\school}{m}{\newcommand{\schoolname}{#1}}
\NewDocumentCommand{\department}{m}{\newcommand{\departmentname}{#1}}
\NewDocumentCommand{\degree}{m}{\newcommand{\degname}{#1}}
\NewDocumentCommand{\course}{m}{\newcommand{\coursename}{#1}}

% Defines my hline ~
\newcommand\myhline{\noindent\rule[.4pt]{\linewidth}{.4pt}\par}

% Glossary & Acronyms
\usepackage[toc,acronym]{glossaries} 
% https://www.dickimaw-books.com/gallery/glossaries-styles/#tree-style
% https://tex.stackexchange.com/questions/346681/space-between-acronym-entry-and-it-description
\setglossarystyle{long} % Glossary style
\setlength\LTleft{-6pt} % Left margin
\setlength\LTright{0pt} % Right margin
\setlength\glsdescwidth{0.94\hsize} % Glossary description width
\renewcommand{\glsnamefont}[1]{\textbf{#1}} % Notation in bold
\renewcommand{\glossaryentrynumbers}[1]{\textit{(p. \textcolor{darkred}{#1})}} % Page number customisation

% Adjust the height of the text area
\addtolength{\textheight}{-20pt}
\setlength{\headheight}{20pt}

% Header & Footer
\usepackage{fancyhdr} % Required for customizing headers and footers
\usepackage[bottom, hang]{footmisc} % Force footnotes to the bottom
\pagestyle{fancy} \fancyhf{} % Clear default style
\renewcommand{\headrulewidth}{.3pt}
\setlength{\footnotemargin}{6pt} % Footnote space

\renewcommand{\chaptermark}[1]{\markboth{\textbf{\thechapter. #1}}{}}
\renewcommand{\sectionmark}[1]{\markright{\textbf{\thesection. #1}}}

\fancyhead[OLH]{
    \ifnum\value{section}=0\leftmark
    \else\rightmark
    \fi
\vspace{3pt}} % Section name... if not empty!
%\fancyhead[L]{\textit{Passer} sp. Genomics}
\fancyhead[L]{Supplementary Information for}
\fancyhead[R]{Pacheco et al. 2025}




% Place the page numbering on bottom when the page style is empty, i.e., when a new chapter starts
\fancypagestyle{plain}{%
  \fancyhf{}%
  \renewcommand{\headrulewidth}{0pt}%
  \fancyfoot[ORF,ELF]{\thepage\vspace{3pt}}%
}

% Bibliography
\usepackage[
	backend=biber, % Use the biber backend for compiling the bibliography
	citestyle=authoryear, % In-text citation style
	bibstyle=authoryear, % Bibliography style
	sorting=nyt, % Order references
    natbib=true, % Use natbib compatibility mode
    dashed=false, % Always print the author in the bibliography
    uniquelist=false % Only show one author, e.g. (Author et al., 2023)
]{biblatex}

\DefineBibliographyStrings{english}{
  bibliography = {Bibliografia (EN)},
}

\renewcommand{\appendixname}{Cartas de Apoio}

\setlength\bibitemsep{1.5\itemsep} % Adjust the space between references
\setlength\bibhang{0pt} % Remove indentation
\renewcommand*{\bibfont}{\small} % Change references font size
\addbibresource{Bibliography/Bibliography.bib} % Reference document

% Highlight both Author and Year in the citations
\ExecuteBibliographyOptions{maxcitenames=1}
\DeclareFieldFormat{citehyperref}{%
  \DeclareFieldAlias{bibhyperref}{noformat}% Avoid nested links
  \bibhyperref{#1}}
\DeclareFieldFormat{textcitehyperref}{%
  \DeclareFieldAlias{bibhyperref}{noformat}% Avoid nested links
  \bibhyperref{%
    #1%
    \ifbool{cbx:parens}
      {\bibcloseparen\global\boolfalse{cbx:parens}}
      {}}}
\savebibmacro{cite}
\savebibmacro{textcite}
\renewbibmacro*{cite}{%
  \printtext[citehyperref]{%
    \restorebibmacro{cite}%
    \usebibmacro{cite}}}
\renewbibmacro*{textcite}{%
  \ifboolexpr{
    ( not test {\iffieldundef{prenote}} and
      test {\ifnumequal{\value{citecount}}{1}} )
    or
    ( not test {\iffieldundef{postnote}} and
      test {\ifnumequal{\value{citecount}}{\value{citetotal}}} )
  }
    {\DeclareFieldAlias{textcitehyperref}{noformat}}
    {}%
  \printtext[textcitehyperref]{%
    \restorebibmacro{textcite}%
    \usebibmacro{textcite}}}
    
% Tables
\usepackage{tabularx} % Table length
\usepackage{booktabs} % Table style
\usepackage{multirow} % Multirow
\usepackage{longtable} % For long tables
\newcommand{\customtableformatting}{
    \renewcommand{\arraystretch}{1.1}\footnotesize
} % Increase space between rows and smaller font
\AtBeginEnvironment{tabular}{\customtableformatting} % Apply to tabular
\AtBeginEnvironment{tabularx}{\customtableformatting} % Apply to tabularx
\preto\longtable{\customtableformatting} % Apply to longtable


\renewcommand{\thefigure}{S\arabic{figure}}

% Captions
\usepackage{caption} % Required for customizing captions
% Set vertical space before every table
\setlength{\abovecaptionskip}{10pt}
\captionsetup{skip=6pt} % Whitespace between figures/tables and the caption
\captionsetup{
    labelfont={bf,small},
    textfont={small}
} % Caption font style
\captionsetup[table]{skip=3pt, justification=centering}
\captionsetup[figure]{justification=justified, name=Supplementary Figure}

% Lists
\usepackage{enumitem} % Required for list customization
\setlist[itemize,enumerate,description]{left=0pt, topsep=0pt, noitemsep}
\setlist{noitemsep} % Customize spacing around and inside lists

% Links
\usepackage[bookmarks]{hyperref} % Required for links

% Capitalise the first letter when \autoref for both English and Portuguese
\addto\extrasenglish{\def\figureautorefname{Figure}}
\addto\extrasenglish{\def\chapterautorefname{Chapter}}
\addto\extrasenglish{\def\sectionautorefname{Section}}
\addto\extrasenglish{\def\subsectionautorefname{Section}}
\addto\extrasenglish{\def\subsubsectionautorefname{Section}}
\addto\extrasenglish{\def\tableautorefname{Table}}
\addto\extrasenglish{\def\partautorefname{Part}}
\addto\extrasenglish{\def\appendixautorefname{Appendix}}
\addto\extrasportuguese{\def\sectionautorefname{Secção}}
\addto\extrasportuguese{\def\subsectionautorefname{Secção}}
\addto\extrasportuguese{\def\subsubsectionautorefname{Secção}}

% Custom color for citations and references
\hypersetup{
	colorlinks=true, % Whether to color the text of links
	urlcolor=[HTML]{045a8d}, % Color for \url and \href links
	linkcolor=darkred, % Color for \nameref links
	citecolor=[HTML]{7f2704}, % Color of reference citations
}

% Chapters & Sections
\usepackage{titlesec}
\titleformat{\chapter}[display]
{\flushright\normalsize\huge\color{black}}%
{\flushright\normalsize%
{\fontsize{25}{25}\selectfont\bfseries\hspace*{-1em}\rule[-2pt]{1.5pt}{30pt}\hspace{0.5em}\raisebox{4.5pt}{\thechapter}}}
{10pt}%
{\addfontfeature{LetterSpace=6}\huge\scshape\bfseries}

\titleformat{\section}
{\large\bfseries} % Format: Large size, bold
{\thesection} % Section number style
{0.75em} % Spacing between number and title
{} % Code before the title (empty in this case)
\titlespacing*{\section} % For \section
{0pt}  % Left margin (indentation)
{12pt} % Space before the section
{8pt}  % Space after the section


\titleformat{\subsection}
{\normalsize\bfseries} % Format: Large size, bold
{\thesubsection} % Section number style
{0.75em} % Spacing between number and title
{} % Code before the title (empty in this case)
\titlespacing*{\subsection} % For \subsection
{0pt}  % Left margin (indentation)
{8pt} % Space before the section
{6pt}  % Space after the section

% Abstract Keywords
\NewDocumentCommand{\keywordsen}{m}{
    \vspace{.4cm}\noindent\textbf{Keywords:} #1
}
\NewDocumentCommand{\keywordspt}{m}{
    \vspace{.4cm}\noindent\textbf{Palavras-Chave:} #1
}

% Blankpage with warning for left blank
\newcommand\blankpage{
    \clearpage
    \thispagestyle{empty}%
    \addtocounter{page}{1}%
    \vspace*{\fill}
    \begin{center}
        \ifthenelse{\equal{\getLanguage}{portuguese}}{%
            \textcolor{gray!50}{\textit{Página intencionalmente deixada em branco.}}
        }{%
            \textcolor{gray!50}{\textit{This page intentionally left blank.}}
        }
    \end{center}
    \vspace*{\fill}
    \clearpage
}

% Defines the colour for the important box background
\definecolor{importantboxcolour}{HTML}{ece7f2 }
  
% Define a custom environment for the highlighted box
\newenvironment{importantbox}{%
  \begin{tcolorbox}[colback=importantboxcolour, colframe=importantboxcolour, arc=5pt, boxrule=1pt, left=15pt, right=15pt, top=10pt, bottom=10pt]%
  \small
}{%
  \end{tcolorbox}%
}

% Defines the colours for title and background
\definecolor{customdarkblue}{HTML}{807dba}
\definecolor{customlightblue}{HTML}{dadaeb}
\definecolor{customdarkgreen}{HTML}{365338}
\definecolor{customlightgreen}{HTML}{7E987D}
\definecolor{customdarkred}{HTML}{ee7576}
\definecolor{customlightred}{HTML}{f0f0f0}


% Define of both \note and \todo commands
\newtcolorbox{custombox}[2][]{enhanced,
    before skip=2mm,after skip=2mm,colframe=gray!30,boxrule=0.2mm,
    attach boxed title to top left={xshift=1cm,yshift*=1mm-\tcboxedtitleheight},
    varwidth boxed title*=-3cm,
    boxed title style={frame code={
        \path[fill=tcbcolback!30!black]
            ([yshift=-1mm,xshift=-1mm]frame.north west)
            arc[start angle=0,end angle=180,radius=1mm]
            ([yshift=-1mm,xshift=1mm]frame.north east)
            arc[start angle=180,end angle=0,radius=1mm];
        \path[left color=tcbcolback!80!black,right color=tcbcolback!80!black,
            middle color=tcbcolback!80!black]
            ([xshift=-2mm]frame.north west) -- ([xshift=2mm]frame.north east)
            [rounded corners=1mm]-- ([xshift=1mm,yshift=-1mm]frame.north east)
            -- (frame.south east) -- (frame.south west)
            -- ([xshift=-1mm,yshift=-1mm]frame.north west)
            [sharp corners]-- cycle;
        },interior engine=empty,
    },
    fonttitle=\MakeUppercase\small\ttfamily\bfseries,
    title={#2},#1
}

\newcommand{\todo}[1]{
    \vspace{10pt}
    \begin{custombox}[colbacktitle=customdarkblue,colback=customlightblue,breakable]{CITATION}
        \small\itshape #1
    \end{custombox}
    \vspace{10pt}
}

\newcommand{\pressing}[1]{
    \vspace{10pt}
    \begin{custombox}[colbacktitle=customdarkgreen,colback=customlightgreen,breakable]{PRESSING QUESTIONS}
        \small\bfseries\textcolor[HTML]{ffffff}{#1}
    \end{custombox}
    \vspace{10pt}
}

\newcommand{\ofnote}[1]{
    \vspace{10pt}
    \begin{custombox}[colbacktitle=customdarkred,colback=customlightred,breakable]{OF NOTE}
        \small\bfseries #1
    \end{custombox}
    \vspace{10pt}
}

\newcommand{\note}[1]{
    \vspace{10pt}
    \begin{custombox}[colbacktitle=blue,colback=blue!10,breakable]{NOTE}
        \small\ttfamily #1
    \end{custombox}
    \vspace{10pt}
}

% Code highlighting and printing
\usepackage[newfloat]{minted}
\newenvironment{code}{\captionsetup{type=listing}}{}
\ifthenelse{\equal{\getLanguage}{portuguese}}{%
    \SetupFloatingEnvironment{listing}{name=Listagem}
}{%
    \SetupFloatingEnvironment{listing}{name=Listing}
}
\captionsetup[listing]{skip=-7pt}
\captionsetup{justification=raggedright,singlelinecheck=false}

\newenvironment{longlisting}{\captionsetup{type=listing}}

% Define a new minted style
\setminted{
    frame=lines,
    framesep=2mm,
    baselinestretch=1.2,
    fontsize=\footnotesize,
    linenos,
    tabsize=4,
    breaklines=true,
    breakanywhere=true
}

% Create a plainpage followed by a blankpage
\newcommand{\plainblankpage}{\thispagestyle{plain}\blankpage}

% Landscape mode done correctly
% https://tex.stackexchange.com/questions/9071/how-to-translate-and-rotate-the-heading-of-landscaped-pages
\makeatletter
\def\ifGm@preamble#1{\@firstofone}
\appto\restoregeometry{%
    \pdfpagewidth=\paperwidth
    \pdfpageheight=\paperheight}
\apptocmd\newgeometry{%
    \pdfpagewidth=\paperwidth
    \pdfpageheight=\paperheight}{}{}
\makeatother

\newenvironment{landscapemode}[2]{%
    \newgeometry{paperwidth=#1,paperheight=#2,hmargin=3cm,vmargin=5cm,top=3cm,landscape}
    \fancyheadoffset{0pt}
    %\titlespacing*{\chapter}{0pt}{0pt}{40pt}
}{%
    \restoregeometry
    \fancyheadoffset{0pt}
}

% Define a new command for "BrainyBirdS"
\newcommand{\BrainyBirds}{\textbf{\textit{\textcolor[HTML]{5d585b}{Brainy}\textcolor[HTML]{827a7f}{Birds}}\textcolor[HTML]{009B77}{\textsuperscript{B}}\textcolor[HTML]{FFDD33}{\textsubscript{R }}}}