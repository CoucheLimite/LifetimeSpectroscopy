# LifetimeSpectroscopy
This is a Python based GUI for lifetime spectroscopy analysis
It can read Sinton excel file and exported txt data file from QSSPL analyzer (https://github.com/MK8J/QSSPL-analyser)
It has some small data pre-process functions like cropping, intrinsic correction etc. but its main job is to do lifetime spectroscopy analysis.
It can analysis data using three different methods:
1. Defect parameter solution surface (DPSS method) using Murphy linearization [1-3]
2. Newton method [4]
3. Simultaneouly fitting [5, 6]
Read the reference below to get more information on these methods

## Installation
Simply download all the files and double click IDLSanalyzer.py or run it with Python

### Prerequisites
1. Install Python 3.6
2. Install numpy, scipy, matplotlib, openpyxl (for windows uuser, it is recommended to get thses packages from http://www.lfd.uci.edu/~gohlke/pythonlibs/ and use pip to install, be careful when choose the version)
3. Install semiconductor from MK8J's Github from https://github.com/MK8J/semiconductor (You can download the whl file in the dist folder using pip)
4. install PyQt5 using pip install PtQt5
5. intall pandas using pip install pandas

It seems the PyQt5 only support python v3.5 and later. Sorry for users of other versions!
The GUI is written and tested in 64-bit platform, but it should work well in 32-bit platform. Please let me know if you met any issues.

## Usage

1. Load the file
2. Set or check parameters. (The color of the file will change to green after setting the parameters. Parameters will be read directly for Sinton excel files.)
3. Process if you want. (merge two files, crop, inverse subtraction, subtract intrinsic)
4. Add interested file to analysis list and start analysis
5. Choose either fitting individually or fitting simultaneouly
6. For fitting individually, choose just one file, select to fit with one defect or two (you can adjust the percentage for fitting the initial value), accept fitting if you are satisfied, choose the fiited defects for DPSS or Newton method 
7. For fitting simultaneouly, choose multiple files and set the fitting range, and click fit. This method can take a long time so start with less steps

Everything should be pretty straightfoward to use. Each pushbuttom is enabled only after all the necessary requirements are met. (correct selection from the list, parameters inputed, ect.)
If you have any problems or suggestions, please feel free to contact me.

## Reference

[1] S. Rein, Lifetime Spectroscopy: A Method of Defect Characterization in Silicon for Photovoltaic Applications (Springer Science and Business Media, Berlin, 2006).

[2] J.D. Murphy, K. Bothe, R. Krain, V. V Voronkov, and R.J. Falster, J. Appl. Phys. 111, 113709 (2012).

[3] A.E. Morishige, M.A. Jensen, D.B. Needleman, K. Nakayashiki, J. Hofstetter, T.A. Li, and T. Buonassisi, IEEE J. Photovoltaics 6, 1466 (2016).

[4] Y. Zhu, Q.T. Le Gia, M.K. Juhl, G. Coletti, and Z. Hameiri, IEEE J. Photovoltaics 1 (2017).

[5] C. Sun, F.E. Rougieux, and D. Macdonald, J. Appl. Phys. 115, 214907 (2014).

[6] S. Bernardini, T.U. Naerland, A.L. Blum, G. Coletti, and M.I. Bertoni, Prog. Photovoltaics Res. Appl. 25, 209 (2017).


## Update log
14/08/2017 Major updates, simutaneously fitting enabled
26/08/2017 Add plot of the standard deviation of the DPSS curves
31/08/2017 MajorEnable to export the data in the plot and the calculation of the standard deviation of the DPSS curves. Fixed some bugs in changing model.

To be completed:
1. <strike> Surface recombination subtraction </strike><br/>
2. <strike> Export data </strike><br/>
3. <strike> Calculation of the standard deviation of the DPSS curves</strike><br/>
4. Complete the coding for advanced fitting mode
5. Add the consideration of measurement uncertainty to calculate Chi-squared error.
6. Add the option to minimize the Chi-squared error. (Currently it minimize the sum of squared error)
7. Add the function to get the acceptance range of fitting result.
