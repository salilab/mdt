\renewenvironment{sidewaystable}{\begin{table}}{\end{table}}

\newcommand{\Options}[1]{\begin{tabular}{lll} 
			#1 
			\end{tabular}}

\newcommand{\OptLine}[5]{
       {{\small\sloppy \K{#1} = #2{#3} }}
     & {{\small\sloppy\tt #4}}
     & {{\small\sloppy #5}}\\ }
