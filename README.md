## tsa4
<a href="https://github.com/nickpoison"><img src="https://img.shields.io/badge/NickyPoison-approved-ff69b4.svg?style=flat"></a>

* Because R code takes up so much space in the text, most of the R code for the *graphics* in Time Series Analysis and Its Applications, Edition 4 was shortened to the minimum number of lines needed to get a decent plot.

* The actual code used in the text *for the graphics* is listed in these xxx.R files.  

* Some of the code was written before `astsa` took full form and was the incentive to write our own scripts.  For example, in Chapter 1, you might see 10 or more lines just to get the ACF and PACF on the same scale... this stuff was written before `acf2` and I never updated  it because it did the job.

* The code for *data analysis* is shown in the text exactly as it was used.  All that stuff is on the [website for the text](http://www.stat.pitt.edu/stoffer/tsa4/).

* `grid.r` as listed here was used instead of writing out the details each time... the change was to `col = gray(.9), lty = 1` in the defaults from the original `col = "lightgray", lty = "dotted"` ... not a big deal but a pain if it's called a lot.

* Chapter 7 code is long, so there are a number of files in that directory.

* The use of `ggplot2` in the new edition of the text was considered for about 2 seconds, but it's not what the text is about and it would have added too much space to address it.  Instead, a web page is devoted to `ggplot`, `ggfortify` and base graphics. See the [graphics fix](http://www.stat.pitt.edu/stoffer/tsa4/tsgraphics.htm) page at [tsa4](http://www.stat.pitt.edu/stoffer/tsa4/).
