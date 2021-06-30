# Coverage_Bias_Detection

Uses a piecewise regression approach to detect bias in coverage for RNA-Seq runs. 

`!python coverage_test_piecewise -i compiled_matrix_5_to_3.txt`

The input is the output for picard CollectRnaSeqMetrics only adding to a separate file the normalized position/normalized coverage and values part of the file.

Program prints out whether the specific sample passed or failed. 
