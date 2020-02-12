# Image Processing

University project for Procedural Programming course. It is an application that encrypts and decrypts a BMP image and does handwritten character recognition using Template Matching.

# Features

- encryption and decryption a BMP image using a symmetrical cipher
- assessment of the uniformity of the distribution of values in a string using the chi-squared test
- recognition of handwritten characters

# Compile and Run

``` 
gcc main.c -o main -lm
./main
```
# Output

### Encryption and Decryption

peppers.bmp           |  peppers_enc.bmp         |   peppers_dec.bmp
:-------------------------:|:-------------------------:|:-------------------------:
<img src = "https://github.com/AMAPuiu/Image-Processing-/blob/master/Images/peppers.bmp" width="200" > | <img src = "https://github.com/AMAPuiu/Image-Processing-/blob/master/Images/peppers_enc.bmp" width="200"> | <img src = "https://github.com/AMAPuiu/Image-Processing-/blob/master/Images/peppers_dec.bmp" width="200">

test.bmp           |  test_rep.bmp        
:-------------------------:|:-------------------------:
<img src = "https://github.com/AMAPuiu/Image-Processing-/blob/master/Images/test.bmp" width="200" > | <img src = "https://github.com/AMAPuiu/Image-Processing-/blob/master/Images/test_rec.bmp" width="200">
