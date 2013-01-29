
%%
%% XXXII Iberian Latin-American Congress on Computational Methods in Engineering (XXXII Cilamce).
%% Ouro Preto, Minas Gerais, Brazil
%% http://www.acquacon.com.br/cilamce2011/en/
%% 
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{cilamce2011}[2011/07/06 Cilamce 2011 paper class]
%
\DeclareOption{10pt}{\OptionNotUsed}
\DeclareOption{11pt}{\OptionNotUsed}
\DeclareOption{twocolumn}{\OptionNotUsed}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
%
\def\cilamcelanguage{1}
\DeclareOption{spanish}{
  \def\cilamcelanguage{0}
}
\def\cilamcepagenumbers{0}
\DeclareOption{pagenumbers}{
  \def\cilamcepagenumbers{1}
}
%
\DeclareOption{english}{
  \def\cilamcelanguage{1}
}
%
\DeclareOption{portuguese}{
  \def\cilamcelanguage{2}
}
%
\def\cilamceuselinks{1}
\DeclareOption{links}{
  \def\cilamceuselinks{1}
}
\DeclareOption{nolinks}{
  \def\cilamceuselinks{0}
}
%
\ProcessOptions\relax
%
\LoadClass[12pt]{article}
%
\RequirePackage{calc}
\RequirePackage{indentfirst}
%
% This allows that copying from PDF files preserves accents. 
\RequirePackage{times}
\RequirePackage[T1]{fontenc} % Usar T1 enconding
%
%% Title
%
\renewcommand{\@maketitle}{%
  \newpage
  \null
  \vspace{2cm}%
  {\begin{center}%
     {\fontsize{14}{16}\selectfont \bfseries \MakeUppercase{\@title} \par}%
     \vspace{12pt}%
     {
        \lineskip .5em%
        {\setlength{\tabcolsep}{0in}
        \begin{tabular}[t]{c}%
          \bfseries \@author
        \end{tabular}\par} }%
  \end{center}}% end of \parbox
  } % end of \@maketitle

%\newcommand{\abstractname}{Abstract}
\renewenvironment{abstract}{%
      \noindent \bfseries \abstractname.%
      \mdseries %\small%
}{\thispagestyle{empty}%
\pagebreak}
%
\newenvironment{keywords}
{\vspace{14pt}\noindent\textbf{\keywordsname: }}
{\vspace{12pt}}
%
\renewenvironment{thebibliography}[1]
     {\section*{\refname
        \@mkboth{\MakeUppercase\refname}{\MakeUppercase\refname}}%
      \list{\@biblabel{\@arabic\c@enumiv}}%
           {\settowidth\labelwidth{\@biblabel{#1}}%
            \leftmargin\labelwidth
            \parsep 0pt
            \itemsep 0pt
            \advance\leftmargin\labelsep
            \@openbib@code
            \usecounter{enumiv}%
            \let\p@enumiv\@empty
            \renewcommand\theenumiv{\@arabic\c@enumiv}}%
      \sloppy
      \clubpenalty4000
      \@clubpenalty \clubpenalty
      \widowpenalty4000%
      \sfcode`\.\@m}
     {\def\@noitemerr
       {\@latex@warning{Empty `thebibliography' environment}}%
      \endlist}
%
\renewcommand{\thesubsection}   {\textup{\thesection.\arabic{subsection}}}
\renewcommand{\thesubsubsection}{\textup{%
    \thesubsection .\arabic{subsubsection}}}
%
\renewcommand{\section}{\@startsection{section}{1}{0mm}%
                                     {-12pt}%
                                     {6pt}%
                                     {\bfseries\normalsize\MakeUppercase}}%
\renewcommand{\subsection}{\@startsection{subsection}{2}{0mm}%
                                     {-12pt}%
                                     {6pt}%
                                     {\bfseries\normalsize}}%

%
\setlength{\textfloatsep}{12pt}
\setlength{\abovedisplayskip}{6pt}
\setlength{\belowdisplayskip}{6pt}
\setlength{\parindent}{5mm}
\setlength{\abovecaptionskip}{6pt}
\setlength{\belowcaptionskip}{6pt}
%
\setlength{\textwidth}{16cm}
\setlength{\textheight}{24cm}
\setlength{\oddsidemargin}{-0.04cm}
\setlength{\topmargin}{8.1mm}
\setlength{\headheight}{0mm}
\setlength{\headsep}{0mm}
%
\renewcommand{\@makecaption}[2]{%
  \vspace{\abovecaptionskip}%
  {\fontsize{10}{12} \selectfont
  \sbox\@tempboxa{#1: #2}%
  \ifdim \wd\@tempboxa >\hsize
    #1: #2\par
  \else
    \global \@minipagefalse
    \hb@xt@\hsize{\hfil\box\@tempboxa\hfil}%
  \fi}
  \vskip\belowcaptionskip}
%
\def\voidaffil{$\mskip 0mu$}
%
\RequirePackage[noblocks]{authblk}
\renewcommand\AB@affilsepx{\protect\vspace{12pt}\\\protect\Affilfont }
\renewcommand\Authands{ and }
\renewcommand{\Affilfont}{\mdseries\itshape\small}
\renewcommand{\Authfont}{\bfseries}
%
\RequirePackage{natbib}
\bibliographystyle{cilamce2011}
%\bibliographystyle{unsrt}
\bibpunct{(}{)}{;}{a}{,}{,}
\setlength{\bibsep}{0pt}
%
\ifnum \cilamcepagenumbers = 0
\pagestyle{empty}
\fi
%
%
\ifnum \cilamcelanguage = 0
  \renewcommand\Authand{ y }%
  \renewcommand\Authands{ y }
  \RequirePackage[spanish]{babel}
  \def\keywordsname{Palabras Clave}
  \def\cilamcebtxand{y}
  \def\cilamcebtxeditors{editores}
  \def\cilamcebtxeditedby{editado por}
  \def\cilamcebtxedition{edici\'on}
  \def\cilamceabtxvolume{volumen}
  \def\cilamcebtxof{de}
  \def\cilamcebtxnumber{n\'umero}
  \def\cilamcebtxpages{p\'aginas}
  \def\cilamcebtxchapter{cap\'\i{}tulo}
  \def\cilamcebtxtechrep{Informe T\'ecnico}
  \def\cilamcebtxmaster{Tesis de Maestr\'\i{}a}
  \def\cilamcebtxphd{Tesis de Doctorado}
  \def\cilamcebtxpage{p\'agina}
  \addto\captionsspanish{%
     \def\tablename{Tabla}}%
\fi
%
\ifnum \cilamcelanguage = 1
  \renewcommand\Authand{ and }%
  \renewcommand\Authands{ and }
  \RequirePackage[english]{babel}
  \def\keywordsname{Keywords}
  \def\cilamcebtxand{and}
  \def\cilamcebtxeditors{editors}
  \def\cilamcebtxeditedby{edited by}
  \def\cilamcebtxedition{edition}
  \def\cilamcebtxvolume{volume}
  \def\cilamcebtxof{of}
  \def\cilamcebtxnumber{number}
  \def\cilamcebtxpages{pages}
  \def\cilamcebtxchapter{chapter}
  \def\cilamcebtxtechrep{Technical Report}
  \def\cilamcebtxmaster{Master's Thesis}
  \def\cilamcebtxphd{D.Sc. Thesis} %\def\cilamcebtxphd{Ph.D. thesis}
  \def\cilamcebtxpage{page}
\fi
%
\ifnum \cilamcelanguage = 2
  \renewcommand\Authand{ e }%
  \renewcommand\Authands{ e }
  \RequirePackage[portuguese]{babel}
  \def\keywordsname{Palavras Chave}
  \def\cilamcebtxand{e}
  \def\cilamcebtxeditors{editores}
  \def\cilamcebtxeditors{editores}
  \def\cilamcebtxeditedby{editado por}
  \def\cilamcebtxedition{edi\c c\~ao}
  \def\cilamcebtxvolume{volume}
  \def\cilamcebtxof{de}
  \def\cilamcebtxnumber{n\'umero}
  \def\cilamcebtxpages{p\'aginas}
  \def\cilamcebtxchapter{cap\'\i{}tulo}
  \def\cilamcebtxtechrep{Relat\'orio T\'ecnico}
  \def\cilamcebtxmaster{Dissertação de Mestrado}
  \def\cilamcebtxphd{Tese de Doutorado}
  \def\cilamcebtxpage{pagina}
\fi
%
\ifnum \cilamceuselinks = 1

% Use color links
\RequirePackage{color}
\definecolor{url}{rgb}{0.8,0.,0.}
\RequirePackage[colorlinks,linkcolor=blue,%
                urlcolor=url,citecolor=blue,%
                plainpages=false,pdfpagelabels]{hyperref}
\RequirePackage{nameref,url}
\def\UrlFont{\mdseries}
\fi

%%End of file