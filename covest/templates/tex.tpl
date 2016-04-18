\documentclass[a4paper,12pt]{article}
\usepackage[utf8]{inputenc}
\usepackage{lmodern}
\usepackage[T1]{fontenc}
\usepackage{fontenc}
\usepackage{booktabs}

\begin{document}

\begin{table}[htp]
\catcode`\-=12
\centering
\begin{tabular}{ {{# header}}c{{/header}} }
\toprule
{{# header }}{{ value }}{{^ last }} & {{/ last }}{{/ header }}\\
\midrule
{{# body }}
{{# line }}{{ value }}{{^ last }} & {{/ last }}{{/ line }}\\
{{/ body }}
\bottomrule
\end{tabular}
\end{table}

\end{document}
