How to compile the program ?
----------------------------

Type
$make 
in command line is enough.


How to run the program ?
------------------------

Copy the exe file named "giscup" into a folder with dataset file, queries file and all data files.

To run the program
$./giscup your_dataset.txt your_queries.txt
or
$./giscup
with default dataset file named "dataset.txt" and default  queries file named "queries.txt".


What the central idea of the program is ?
-----------------------------------------

This program mainly works on finding similar trajectories in a short time using a mathematics concept - Fréchet Distance. It implements the algorithm to decide whether the Fréchet Distance of two given trajectories is not greater than the given bound. R-tree is used as a spatial access method and parallel computing is also included in this program.