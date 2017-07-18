/*
Copyright 2017 Paul Reggentin

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


*/

math = require('./kmath.js');
//require("./complex.js");
fs = require("fs");
require('./nifti.js');


var filename = 'rotated.nii'

var raw = fs.readFileSync(filename);

var nii = parse(raw);

var orientation = new Float32Array(3);

orient(nii.data,nii.sizes,orientation);

console.log(orientation);
/*
returns information about the orientation of an object in an image
assumes only one object is present
ARGUMENTS:
            img: a float vector representing an image stack
            dim: an integer array representing the dimensions of img
            orientation: a 3-element array that will store the output information
                         [center_x, center_x, angle_of_orientation]
*/
function orient(img, dim, orientation){
    var st = performance.now();

    var slicelen = dim[0]*dim[1];
//binarize image
    //find min and max of histogram by finding extremes of a middle slice
    var middle_slice = Math.floor(dim[2] / 2);
    var min = img[middle_slice], max = img[middle_slice];
    for(var i = middle_slice*slicelen; i < (middle_slice+1)*slicelen; i++){
        if(img[i] < min) {
            min = img[i];
        } else if(img[i] > max){
            max = img[i];
        }
    }
    var range = max - min;
    max += .1*range;
    min -= .1*range;
    var bin_width = (max - min)/254;
    var hist = new Float32Array(256);
    //for each pixel in img, find what bin it belongs to and increment that one
    for(var i = 0; i < img.length; i++){
        var idx = Math.floor((img[i]-min)/bin_width);
        if(idx<0){ idx = 0; }
        else if (idx > 255) { idx = 255; }
        hist[idx] += 1;
    }

    var thresh = min + bin_width * otsu(hist,img.length);
    var binary_img = new Uint8Array(img.length);
    for(var col = 0; col < dim[1]*dim[2]; col++){
        var firstpix = 0;
        var lastpix = dim[0]-1;
        while(firstpix<lastpix && img[dim[0]*col+firstpix] < thresh){
            firstpix++;
        }
        if(firstpix==lastpix){ continue; }
        while(firstpix<lastpix && img[dim[0]*col+lastpix] < thresh){
            firstpix--;
        }
        for(var i = firstpix; i <= lastpix; i++){
            binary_img[dim[0]*col+i] = 1;
        }
    }

    var en = performance.now();
    console.log(en-st);
    exportimg(binary_img,'/home/paul/Documents/inertia_tensor/bin.csv');
    
//find the center of mass

//find the orientation via the principal inertial axes
}

//total is total number of pixels
function otsu(histogram, total) {
    var sum = 0;
    for (var i = 0; i < histogram.length; ++i) //normally it will be 255 but sometimes we want to change step
        sum += i * histogram[i];
    var sumB = 0;
    var wB = 0;
    var wF = 0;
    var mB;
    var mF;
    var max = 0.0;
    var between = 0.0;
    var threshold1 = 0.0;
    var threshold2 = 0.0;
    for (var i = 0; i < histogram.length; ++i) {
        wB += histogram[i];
        if (wB == 0)
            continue;
        wF = total - wB;
        if (wF == 0)
            break;
        sumB += i * histogram[i];
        mB = sumB / wB;
        mF = (sum - sumB) / wF;
        between = wB * wF * (mB - mF) * (mB - mF);
        if ( between >= max ) {
            threshold1 = i;
            if ( between > max ) {
                threshold2 = i;
            }
            max = between;            
        }
    }
    return ( threshold1 + threshold2 ) / 2.0;
}

function exportimg(vec,filename){
    fs.writeFileSync(filename,vec.toString());
}


