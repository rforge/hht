#/bin/bash
R CMD build hht
mv hht_1.0.0.tar.gz tests
cd tests
tar -xvf hht_1.0.0.tar.gz
