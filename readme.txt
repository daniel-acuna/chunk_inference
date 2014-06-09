This package implements the method described in "Multi-faceted aspects of 
chunking enable robust algorithms" by Daniel E. Acuna, Nicholas F. Wymbs,
Chelsea A. Reynolds, Nathalie Picard, Robert S. Turner, Peter L. Strick,
Scott Grafton, and Konrad Kording, under review in the Journal of
Neurophysiology.

The algorithm was implemented by Daniel E. Acuna (daniel.acuna@northwestern.edu)
If you have any questions, send me an email.

To run the algorithm, you need to first detrend you data so that you remove
aspects of training that are largely unrelated to chunking. In our paper,
we described a simple model that does this detrending.

We provide a full worked out example using data from one subject and
one sequence. Our code needs the Statistics Toolbox from Matlab for
detrending.

