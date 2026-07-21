{smcl}
{* *! version 1.0.0 24dec2014}{...}
{hi:help verswlr}{right: ({browse "http://www.stata-journal.com/article.html?article=st0449":SJ16-3: st0449})}
{hline}

{title:Title}

{p2colset 5 16 18 2}{...}
{p2col :{hi:verswlr} {hline 2}}Perform versatile weighted log-rank test{p_end}
{p2colreset}{...}


{title:Syntax}

{p 8 17 2}
{cmd:verswlr} {varname} {ifin} [{cmd:,}
{it:{help verswlr##options:options}}] 

{p 4 6 2}
{varname} should contain the group indicator variable and must be numeric.

{p 4 6 2}
{cmd:by} is allowed; see {manhelp by D}.{p_end}

{p 4 6 2}
You must {cmd:stset} your data before using {cmd:verswlr}; see
{manhelp stset ST}.


{title:Description}

{pstd}{cmd:verswlr} determines the significance level of the maximum of
G(0,0), G(1,0), and G(0,1) weighted log-rank statistics or three other
user-defined tests for the comparison of two survival curves.  Output and
returned scalars are the sample size, chi-squared statistic for each test,
maximum Z statistic (absolute value), and two-sided p-value.


{marker options}{...}
{title:Options}

{phang}
{opt rho1(#)} and {opt gamma1(#)} specify the weights for the first test.

{phang}
{opt rho2(#)} and {opt gamma2(#)} specify the weights for the second test.

{phang}
{opt rho3(#)} and {opt gamma3(#)} specify the weights for the third test.

{pstd}
The default values are (0 0), (1 0), and (0 1), respectively.


{title:Examples}

{phang}{stata "use stablein, clear" : {bf:. use stablein}}{p_end}
{phang}{stata "stset stime, failure(indicos)" : {bf:. stset stime, failure(indicos)}}{p_end}
{phang}{stata "verswlr trt" : {bf:. verswlr trt}}{p_end}
{phang}{stata "verswlr trt, rho2(2) gamma3(2)"  : {bf:. verswlr trt, rho2(2) gamma3(2)}}{p_end}


{title:Stored results}

{pstd}{cmd:verswlr} stores the following in {hi:r()}:

{synoptset 16 tabbed}{...}
{p2col 5 16 20 2: Scalars}{p_end}
{synopt:{cmd:r(sampsize)}}sample size {p_end}
{synopt:{cmd:r(fh11)}}G(rho1,gamma1) chi-squared statistic{p_end}
{synopt:{cmd:r(fh22)}}G(rho2,gamma2) chi-squared statistic{p_end}
{synopt:{cmd:r(fh33)}}G(rho3,gamma3) chi-squared statistic{p_end}
{synopt:{cmd:r(maxz)}}maximum Z statistic{p_end}
{synopt:{cmd:r(pval)}}two-sided p-value{p_end}


{title:Author}

{pstd}Theodore G. Karrison, University of Chicago{break}
tkarrison@health.bsd.uchicago.edu{p_end}


{title:Also see}

{p 4 14 2}
Article:  {it:Stata Journal}, volume 16, number 3: {browse "http://www.stata-journal.com/article.html?article=st0449":st0449}
{p_end}
