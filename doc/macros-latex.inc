% Calpha, Cbeta, Cgamma:
\newcommand{\Calpha}{$\mbox{C}_\alpha$}
\newcommand{\Cbeta}{$\mbox{C}_\beta$}
\newcommand{\Cgamma}{$\mbox{C}_\gamma$}
\newcommand{\Cgammaone}{$\mbox{C}_{\gamma 1}$}
\newcommand{\Cgammatwo}{$\mbox{C}_{\gamma 2}$}
\newcommand{\Cdelta}{$\mbox{C}_\delta$}
\newcommand{\Cdeltaone}{$\mbox{C}_{\delta 1}$}
\newcommand{\Cdeltatwo}{$\mbox{C}_{\delta 2}$}
\newcommand{\Cepsilonone}{$\mbox{C}_{\epsilon 1}$}
\newcommand{\Cepsilontwo}{$\mbox{C}_{\epsilon 2}$}
\newcommand{\Czeta}{$\mbox{C}_\zeta$}
\newcommand{\Obeta}{$\mbox{O}_\beta$}
\newcommand{\Ogamma}{$\mbox{O}_\gamma$}
\newcommand{\Ogammaone}{$\mbox{O}_{\gamma 1}$}

% Subscript, superscript
% the inclusion of \scriptsize command crashes LaTeX when the macro is used
% in caption and some other environments; this can be corrected by the
% \protect command; The argument of the \subs or \sups command is typeset 
% in text-mode (and can only be used in text-mode (add another \mbox around
% the command to use it in both text and math modes, but that will 
% crash even more). Scriptsize command is relative to the \normalsize
% only. So, if you want to get proper size subscripts with other
% size basic text, include the size for scripts into the argument.
\newcommand{\subs}[1]{${}_{\mbox{\protect\scriptsize #1}}$}
\newcommand{\sups}[1]{${}^{\mbox{\protect\scriptsize #1}}$}

% Degrees Celsius, degrees, Subsite, Subsite', Position, Position'
% In these commands, the argument is processed in math-mode and the
% size of super/sub-scripts will be fine whatever the size of the
% basic font; note the use of the mathmode only style \scriptstyle here,
% as opposed the \scriptsize used in text-mode; explicit \scriptstyle
% in math-mode prevents the error message 'scriptscriptfont undefined'.
\newcommand{\oC}{\mbox{${}^{\scriptstyle\circ}$C}}
\newcommand{\degr}{\mbox{${}^{\scriptstyle\circ}$}}

\newcommand{\ASGL} {\htmladdnormallink{\mbox{\sc Asgl}}
           {http://guitar.rockefeller.edu/asgl/asgl.html}}

% Latin words:
\newcommand{\Denovo}{\mbox{\em De novo}}
\newcommand{\denovo}{\mbox{\em de novo}}
\newcommand{\insitu}{\mbox{\em in situ}}
\newcommand{\Insitu}{\mbox{\em In situ}}
\newcommand{\Abinitio}{\mbox{\em Ab initio}}
\newcommand{\abinitio}{\mbox{\em ab initio}}
\newcommand{\versus}{\mbox{\em versus}}
\newcommand{\etal}{\mbox{\em et~al.}}
\newcommand{\etc } {\mbox{\em etc}}
\newcommand{\cf } {\mbox{\em cf}}
\newcommand{\ie } {\mbox{\em i.e.}}
\newcommand{\eg } {\mbox{\em e.g.}}
\newcommand{\via } {\mbox{\em via}}
\newcommand{\adhoc } {\mbox{\em ad~hoc}}
\newcommand{\invivo } {\mbox{\em in~vivo}}
\newcommand{\invacuo} {\mbox{\em in~vacuo}}
\newcommand{\invitro } {\mbox{\em in~vitro}}
\newcommand{\viceversa } {\mbox{\em vice~versa}}
\newcommand{\aposteriori } {\mbox{\em a~posteriori}}
\newcommand{\apriori } {\mbox{\em a~priori}}
\newcommand{\cis } {\mbox{\em cis}}
\newcommand{\trans } {\mbox{\em trans}}
\newcommand{\dd}{ \; \mbox{d}}
