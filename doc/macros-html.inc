% Calpha, Cbeta, Cgamma:
\newcommand{\Calpha}{${C}_\alpha$}
\newcommand{\Cbeta}{${C}_\beta$}
\newcommand{\Cgamma}{${C}_\gamma$}
\newcommand{\Cgammaone}{${C}_{\gamma 1}$}
\newcommand{\Cgammatwo}{${C}_{\gamma 2}$}
\newcommand{\Cdelta}{${C}_\delta$}
\newcommand{\Cdeltaone}{${C}_{\delta 1}$}
\newcommand{\Cdeltatwo}{${C}_{\delta 2}$}
\newcommand{\Cepsilonone}{${C}_{\epsilon 1}$}
\newcommand{\Cepsilontwo}{${C}_{\epsilon 2}$}
\newcommand{\Czeta}{${C}_\zeta$}
\newcommand{\Obeta}{${O}_\beta$}
\newcommand{\Ogamma}{${O}_\gamma$}
\newcommand{\Ogammaone}{${O}_{\gamma 1}$}

% Subscript, superscript
% the inclusion of \scriptsize command crashes LaTeX when the macro is used
% in caption and some other environments; this can be corrected by the
% \protect command; The argument of the \subs or \sups command is typeset 
% in text-mode (and can only be used in text-mode (add another \mbox around
% the command to use it in both text and math modes, but that will 
% crash even more). Scriptsize command is relative to the \normalsize
% only. So, if you want to get proper size subscripts with other
% size basic text, include the size for scripts into the argument.
 \newcommand{\subs}[1]{${}_{#1}$}
 \newcommand{\sups}[1]{${}^{#1}$}

% Degrees Celsius, degrees, Subsite, Subsite', Position, Position'
% In these commands, the argument is processed in math-mode and the
% size of super/sub-scripts will be fine whatever the size of the
% basic font; note the use of the mathmode only style \scriptstyle here,
% as opposed the \scriptsize used in text-mode; explicit \scriptstyle
% in math-mode prevents the error message 'scriptscriptfont undefined'.
\newcommand{\oC}{${}^{o}$C}
\newcommand{\degr}{${}^{o}$}

\newcommand{\ASGL} {\htmladdnormallink{\mbox{\sc Asgl}}
           {http://guitar.rockefeller.edu/asgl/asgl.html}}

% Latin words:
\newcommand{\Denovo}{{\em De novo}}
\newcommand{\denovo}{{\em de novo}}
\newcommand{\insitu}{{\em in situ}}
\newcommand{\Insitu}{{\em In situ}}
\newcommand{\Abinitio}{{\em Ab initio}}
\newcommand{\abinitio}{{\em ab initio}}
\newcommand{\versus}{{\em versus}}
\newcommand{\etal}{{\em et~al.}}
\newcommand{\etc } {{\em etc}}
\newcommand{\cf } {{\em cf}}
\newcommand{\ie } {{\em i.e.}}
\newcommand{\eg } {{\em e.g.}}
\newcommand{\via } {{\em via}}
\newcommand{\adhoc } {{\em ad~hoc}}
\newcommand{\invivo } {{\em in~vivo}}
\newcommand{\invacuo} {{\em in~vacuo}}
\newcommand{\invitro } {{\em in~vitro}}
\newcommand{\viceversa } {{\em vice~versa}}
\newcommand{\aposteriori } {{\em a~posteriori}}
\newcommand{\apriori } {{\em a~priori}}
\newcommand{\cis } {{\em cis}}
\newcommand{\trans } {{\em trans}}
\newcommand{\dd}{ \; {d}}
