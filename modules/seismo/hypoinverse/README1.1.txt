Notes for Hypoinverse 2000, version 1.1, March 2, 2007

This README1.1 file is a supplement to the original README1.0 file. There
is now a revised release in which the program, source code,
documentation, and new test files are all self-consistent and up to
date.

Version 1.1 adds the ability to fix the origin time (see the terminator
card format), initial distance weighting for the first iterations (DI1
command), label earthquakes with the domain and version of processing
used (VER command), labeling phases with an alternate 3-letter
component code (copied from the station card to the archive line),
translate CUSP digitizer codes to data source codes (DID and DIG
commands), and some complexities for assigning and labeling amplitude
magnitudes as either MX (velocity) or ML (local) (see the XMT command
and new fields at the end of the station archive format lines). 

Version 1.1 also outputs "X" characters in cols 119 and 120 of the
archive station lines if the duration or amplitude magnitudes
(respectively) were not used in the event magnitude. This feature is
useful because the final magnitude weight is not output as a number.
The final magnitude is product of the weight on the phase line and the
weight factors from the station card and the time-dependent weight-out
features from the magnitude correction files.

Version 1.1 compiles successfully under the g77 compiler (as well as
the sun f77 compiler), which I understand means the code is linux
compatible.

A new version (1.1) of the documentation is available as a Word doc
file and as a pdf file. Changes include the new commands and the
additional fields in the output formats.

There are now tar files available for the Hypoinverse source code and
data test files. You can also examine and download individual files
from the website. All files are available at

ftp://ehzftp.wr.usgs.gov/klein/hyp2000

A complete download consists of these files

README1.0		Notes from the v. 1.0 2002-2005 release, revised.
README1.1		This file from 2/2007.
hypyymmdd.tar		Source files, including general purpose subroutines.
			tar filename has year, month & day of its creation.
hyptest.tar		3 versions of test data files, including test3,
			an example of the current 2/2007 USGS-NCSN processing.
	Complete user guide in two different formats:
ftp://ehzftp.wr.usgs.gov/klein/hyp2000/docs/hyp2000-1.1.doc   Word file
ftp://ehzftp.wr.usgs.gov/klein/hyp2000/docs/hyp2000-1.1.pdf   pdf file


