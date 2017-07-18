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
require('./Matrix.js');
numeric = require('./numeric-1.2.6.js')

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
    var mask = new Uint8Array(img.length);
    for(var col = 0; col < dim[1]*dim[2]; col++){
        var firstpix = 0;
        var lastpix = dim[0]-1;
        while(firstpix<lastpix && img[dim[0]*col+firstpix] < thresh){
            firstpix++;
        }
        if(firstpix==lastpix){ continue; }
	
        while(firstpix<lastpix && img[dim[0]*col+lastpix] < thresh){
            lastpix--;
        }
	
        for(var i = firstpix; i <= lastpix; i++){
            mask[dim[0]*col+i] = 1;
        }
    }

    
    //exportimg(binary_img,'/home/paul/Documents/inertia_tensor/bin.csv');
    var tot_mass = 0;
    var cntr_mass = new Float64Array(3);
//find the center of mass
    for(var c = 0; c < dim[2]; c++){
        for(var b = 0; b < dim[1]; b++){
            for(var a = 0; a < dim[0]; a++){
                var idx = [get1Dindex(dim,a,b,c)];
                if(mask[idx]){
                    tot_mass += img[idx];
                    cntr_mass[0] += img[idx]*a;
                    cntr_mass[1] += img[idx]*b;
                    cntr_mass[2] += img[idx]*c;
                }
            }
        }
    }
    cntr_mass[0] /=tot_mass;
    cntr_mass[1] /=tot_mass;
    cntr_mass[2] /=tot_mass;

    
//find the orientation via the principal inertial axes
    I = [[0,0,0],[0,0,0],[0,0,0]];
    for(var c = 0; c < dim[2]; c++){
        for(var b = 0; b < dim[1]; b++){
            for(var a = 0; a < dim[0]; a++){
                var idx = [get1Dindex(dim,a,b,c)];
                if(mask[idx]){
                    var aa = a-cntr_mass[0];
                    var bb = b-cntr_mass[1];
                    var cc = c-cntr_mass[2];
                    I[0][0] += img[idx] * (bb*bb + cc*cc);
                    I[0][1] += img[idx] * (-aa) * bb;
                    I[0][2] += img[idx] * (-aa) * cc;
                    I[1][0] += img[idx] * (-aa) * bb;
                    I[1][1] += img[idx] * (aa*aa + cc*cc);
                    I[1][2] += img[idx] * (-bb) * cc;
                    I[2][0] += img[idx] * (-aa) * cc;
                    I[2][1] += img[idx] * (-bb) * cc;
                    I[2][2] += img[idx] * (aa*aa + bb*bb);
                }
            }
        }
    }
    //DEBUG: normalize I
    /*
    var Itot = 0;
    for(var i = 0; i < 3; i++){
        for(var j = 0; j < 3; j++){
            Itot += I[i][j];
    }
    }

    for(var i = 0; i < 3; i++){
        for(var j = 0; j < 3; j++){
            I[i][j] /= Itot;
    }
    }
    */
    ev = numeric.eig(I);

    //smallest eigenvalue corresponds to primary axis of rotation, which should be anterior-posterior axis
    var eigval = ev.lambda.x;
    var axis;
    if(eigval[0] < eigval[1]){
        if(eigval[0] < eigval[2]){
            axis = 0;
        } else {
            axis = 2;
        }
    } else {
        if(eigval[1] < eigval[2]){
            axis = 1;
        } else {
            axis = 2;           
        }
    }


    var en = performance.now();
    console.log(en-st);


    var ex = img.slice(35*slicelen,36*slicelen);
    var scale = ev.E.x[0][axis]/ev.E.x[1][axis];
    for(var x = 0; x < dim[1]; x++){
        var y = math.round((x-cntr_mass[1])*scale+cntr_mass[0]);
        ex[y+dim[0]*x] = 2*max;
    }
    exportimg(ex,'/home/paul/Documents/inertia_tensor/c.csv');

    orientation[0] = cntr_mass[0];
    orientation[1] = cntr_mass[1];
    orientation[2] = Math.atan(scale);
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


function get1Dindex(dim,a,b,c){
    return (a+dim[0]*b+c*dim[0]*dim[1]);
}

/*
function testspeed(){
    var prod;
    var s = performance.now();
    for(var i = 0; i < 10000000000; i++){
        prod = -i*3.56;
    }
    var e = performance.now();
    console.log('-a*b',e-s);

    var s = performance.now();
    for(var i = 0; i < 10000000000; i++){
        prod = (0-i)*3.56;
    }
    var e = performance.now();
    console.log('(0-a)*b',e-s);

    var s = performance.now();
    for(var i = 0; i < 10000000000; i++){
        prod = -1*i*3.56;
    }
    var e = performance.now();
    console.log('-1*a*b',e-s);
}

testspeed();
*/