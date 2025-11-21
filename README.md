\documentclass[a4paper,11pt]{article}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{geometry}
\usepackage{listings}
\usepackage{xcolor}
\usepackage{hyperref}
\usepackage{booktabs}
\usepackage{amsmath}

% Page Setup
\geometry{margin=1in}
\setlength{\parindent}{0pt}
\setlength{\parskip}{0.5em}

% Code Style
\definecolor{codegray}{gray}{0.9}
\lstset{
    backgroundcolor=\color{codegray},
    basicstyle=\ttfamily\small,
    breaklines=true,
    frame=single,
    framerule=0pt,
    framesep=5pt,
    columns=fullflexible
}

\title{\textbf{XH-pi}: XH-pi interactions detector}
\author{}
\date{}

\begin{document}

\maketitle

\section*{Introduction}
\textbf{XH-pi} is a high-performance CLI tool designed for the automated detection of XH--$\pi$ interactions in protein structures (PDB/CIF). It features an integrated topology preparation pipeline (using Gemmi), supports recursive directory scanning, and utilizes efficient geometric algorithms based on \textit{Hudson} and \textit{Plevin} criteria.

\section*{Features}
\begin{itemize}
    \item \textbf{Automated Hydrogen Addition}: In-memory topology preparation via Gemmi.
    \item \textbf{High Performance}: Parallel processing with multi-core support.
    \item \textbf{Flexible Input}: Handles single files or recursive directory scans.
    \item \textbf{Dual Criteria}: Supports both Hudson and Plevin geometric definitions.
    \item \textbf{Industrial Output}: JSON (default) or CSV formats; customizable output directories.
\end{itemize}

\section*{Installation}
Ensure you have Python 3.9+ installed.

\begin{lstlisting}[language=bash]
git clone https://github.com/yourusername/xhpi.git
cd xhpi
pip install .
\end{lstlisting}

\section*{Quick Start}

\subsection*{1. Configure Environment (Optional)}
Set the default path to your Monomer Library (e.g., CCP4 monomers) for permanent usage.
\begin{lstlisting}[language=bash]
xhpi --set-mon-lib /path/to/monomers/
\end{lstlisting}

\subsection*{2. Run Analysis}
Process all \texttt{.cif} files in a directory recursively and save results to JSON.
\begin{lstlisting}[language=bash]
xhpi ./data_directory
\end{lstlisting}

\subsection*{3. Advanced Usage}
Process specific files, use 8 cores, apply "Remove & Re-add" hydrogen mode, and export as CSV.
\begin{lstlisting}[language=bash]
xhpi ./data/*.cif --jobs 8 --h-mode 3 --file-type csv --out-dir ./results
\end{lstlisting}

\section*{Parameters}

\begin{table}[h]
\centering
\begin{tabular}{@{}ll@{}}
\toprule
\textbf{Option} & \textbf{Description} \\ \midrule
\texttt{<PATH>} & Input files (.cif) or directories to scan recursively. \\
\texttt{--set-mon-lib <DIR>} & Set the default Monomer Library path permanently. \\
\texttt{--mon-lib <DIR>} & Override the library path for the current run. \\
\texttt{--out-dir <DIR>} & Specify a central directory for output files. \\
\texttt{--file-type <STR>} & Output format: \texttt{json} (default) or \texttt{csv}. \\
\texttt{--jobs <NUM>} & Number of parallel CPU cores (Default: 1). \\
\texttt{--h-mode <NUM>} & Hydrogen handling strategy (see below). \\ \bottomrule
\end{tabular}
\end{table}

\subsection*{Hydrogen Modes (--h-mode)}
Default is \textbf{4}.
\begin{itemize}
    \item \texttt{0}: NoChange (Keep existing H)
    \item \texttt{1}: Shift (Move H to standard positions)
    \item \texttt{2}: Remove (Remove all H)
    \item \texttt{3}: ReAdd (Remove and re-add all)
    \item \texttt{4}: ReAddButWater (Re-add all except water)
    \item \texttt{5}: ReAddKnown (Only for known residues)
\end{itemize}

\section*{License}
MIT License.

\end{document}