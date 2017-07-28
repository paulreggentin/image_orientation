/*
Copyright 2017 Paul Reggentin

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/

var numeric = require('./numeric-1.2.6.js')
var FFT = require('fft');


var img_proc = {}

module.exports = img_proc;

/*
Javascript Image Processing Toolbox

FUNCTIONS:
downsample  (DONE)
find_symmetry_plane (DONE)
find_center_hough
otsu        (DONE)
exportimg   (DONE)
fft_2d      (DONE)
ifft_2d     (DONE)
smooth      (IN PROGRESS)

OBJECTS:
CoordTransform (DONE)
IndexTransform (DONE)
N_Dim_Gauss    (DONE)
*/


/*
downsample an image or image stack by a specified factor
arguments:
            old_im: image to be downsampled. stored as 1-D vector
            dim:    dimensions of old_im
            factor: downsampling factor. if scalar, all downsampling factors are uniform
                    if vector, must be same length as dim, and scales each dimension separately
            new_im: vector to store downsampled image in. it's a 1-D vector
*/
img_proc.downsample = function(old_im,dim,factor,new_im){
    //if the factor argument is a scalar, convert it to a vector
    var tfac;

    if(factor.length == undefined){
        tfac = new Int16Array(dim.length);
        for(var i = 0; i < tfac.length; i++){
            tfac[i] = factor;
        }
    } else {
        tfac = factor;
    }

    //make sure the scaling is possible
    for(var i = 0; i < dim.length; i++){
        if(dim[i] % tfac[i] != 0){
            throw 'Error: downsampling factor does not divide image dimension';   
        }
    }

    if(dim.length == 2){
        downsample_single_slice(old_im,0,[dim[0],dim[1],1], [tfac[0], tfac[1], 1],new_im, 0);
    } else {

        old_slicelen = dim[0]*dim[1];
        new_slicelen = old_slicelen / (tfac[0]*tfac[1]);
        for(var i = 0; i < dim[2]/tfac[2]; i++){
            downsample_single_slice(old_im, i*old_slicelen*factor[2], dim, tfac, new_im, i*new_slicelen);
        }
    }

    var scale = 1/(factor[0]*factor[1]*factor[2]);
    for(var i = 0; i < new_im.length; i++){
        new_im[i] *= scale;
    }
}

function downsample_single_slice(old_im,old_offset,dim,factor,new_im, new_offset){
    var big_idx = new img_proc.IndexTransform(dim);
    var lil_idx = new img_proc.IndexTransform([dim[0]/factor[0],dim[1]/factor[1],dim[2]/factor[2]]);

    

    //for each pixel in old_im
    for(var i = 0; i < dim[0]*dim[1]*factor[2]; i++){
        //add its value to the corresponding pixel of new_im
        var c1  = big_idx.get_coord(i+old_offset);
        for(var k = 0; k < 3; k++){
            c1[k] = Math.floor(c1[k]/factor[k]);
        }
        var new_i = lil_idx.get_index(c1);
        if(old_im[i]!=0){
            var q = new_i;
        }

        new_im[new_i] += old_im[i+old_offset];
    }
        
    
    //divide each element of new_im by the number of pixels added to it
    /*
    var scale = 1/(factor[0]*factor[1]*factor[2]);
    var newlen = dim[0]*dim[1]*dim[2]/(factor[0]*factor[1]*factor[2]);
    for(var i = 0; i < newlen; i++){
        new_im[new_offset+i] *= scale;
    }
    */
    

    /*
    for(var c = 0; c<factor[2]; c++){
        for(var b = 0; b<factor[1]; b++){
            for(var a = 0; a<factor[0]; a++){
                
            }
        }
    }
    */
}


//k is the size of the square gaussian filter
img_proc.smooth = function(img,dim,filter_dim){
    //create a 3-D gaussian filter with the dimensions factor
    var filter_idx = new img_proc.IndexTransform(factor);
    var filter = new Float32Array(factor[0]*factor[1]*factor[2]);
    var filter_center = [(factor[0]-1)/2,(factor[1]-1)/2,(factor[2]-1)/2];
    var filter_func = new img_proc.N_Dim_Gauss(filter_center,factor);
    
    var weight = 0;
    for (var i = 0; i < filter.length; i++){
        filter[i] = filter_func.val(filter_idx.get_coord(i));
        weight += filter[i];
    }

    for (var i = 0; i < filter.length; i++){
        filter[i] /= weight;
    }
    img_proc.exportimg(filter,'filter.csv');
}

/*
given three axes, determines which one is the anterior-posterior axis
ARGU;MENTS:
        img: the image, stored as a 1-d vector
        mask: a binary matrix, mask[i] == true if img[i] is part of the brain tissue
        img_to_rw: a CoordTransform object that takes voxel coordinate and transforms it to real world
        index: an IndexTransform object that converts between img index and img coordinate
        axes: the three pricipal axes of the brain (usually found with orient.js)
        center: the center of mass of the brain, also probably found with orient.js
*/
img_proc.find_symmetry_plane = function (img,mask,img_to_rw,index,axes,center){
    //for each axis
    var max_refl = 0;
    var best_plane = -1;
    for(var ax = 0; ax < 3; ax++)
    {
        var this_refl = 0;
        var plane = new img_proc.Plane(axes[ax],center)
        //for each point in mask
        for(var i = 0; i < img.length; i++){
            if(mask[i]){
                var r = plane.reflect(img_to_rw.trans(index.get_coord(i)));
                for(var d = 0; d < 3; d++){
                    r[d] = Math.round(r[d]);
                }
                var r_idx = index.get_index(r);

                var refl_val = img[r_idx];

                if(isNaN(refl_val) ){
                   refl_val = 0;
                }
                this_refl += img[i]*refl_val;
            }
        }

        if(this_refl > max_refl){
           max_refl = this_refl;
           best_plane = ax;
        }
        //do some sort of scaling
    }
    return best_plane;
}


/*
finds the center of a circular or ellipse
*/
img_proc.find_center_hough = function (img,dim,center){
    var index = new img_proc.IndexTransform(dim);

}








//total is total number of pixels
img_proc.otsu = function (histogram, total) {
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

img_proc.exportimg = function (vec,filename){
    fs.writeFileSync('/home/paul/Documents/inertia_tensor/'+filename,vec.toString());
}





/*
    perform a 2-dimensional fourier transform 
    img: the image to be transformed, stored as a 1D vector. 
        can be real or complex, stored as [real0 imag0 real1 imag1 ...]
    dim: the dimensions of the 2D image that 'img' represents
    img_ft: the vector that will hold the fourier transform of img
    is_cplx: boolean flag for type of img
*/
img_proc.fft_2d = function (img, dim, img_ft, is_cplx){
    //create fft objects
    fft0 = new FFT.complex(dim[0],false);
    fft1 = new FFT.complex(dim[1],false);

    //create vectors that will store intermediate steps
    var frow = new Float64Array(2*dim[0]);
    var intermed = new Float64Array(2*dim[0]*dim[1]);
    var intermed_t = new Float64Array(2*dim[0]*dim[1]);
    var fft_trans = new Float64Array(2*dim[0]*dim[1]);

    //assign the correct input row length
    var fft_type;
    var row_len;
    if(is_cplx){
        fft_type = 'complex';
        row_len = 2*dim[0];
    } else {
        fft_type = 'real';
        row_len = dim[0];
    }

    //perform 1D fft on each of the rows of img
    for(var i = 0; i < dim[1]; i++){
        fft0.simple(frow,img.slice(i*row_len,(i+1)*row_len),fft_type);
        copy_into(frow,intermed,2*dim[0],i*2*dim[0]);
    }

    //transpose the row transformed array
    transpose_cplx(intermed,dim,intermed_t);

    //fft each of the rows (which are columns from the original image)
    var col_len = 2*dim[1];
    var fcol = new Float64Array(col_len);
    for(var i = 0; i < dim[0]; i++){
        fft1.simple(fcol,intermed_t.slice(i*col_len,(i+1)*col_len),'complex');
        copy_into(fcol,fft_trans,col_len,i*col_len);
    }

    //transpose back to original dimentions
    dim_t = new Float64Array([dim[1], dim[0]]);
    transpose_cplx(fft_trans,dim_t,img_ft);
}

//inverse fourier transform
//same as fft_2d but all of the 1D ffts are inverse DFTs
img_proc.ifft_2d = function (kspace, dim, img_space){
    fft0 = new FFT.complex(dim[0],true);
    fft1 = new FFT.complex(dim[1],true);

    //fft each of the rows
    var frow = new Float64Array(2*dim[0]);
    var intermed = new Float64Array(2*dim[0]*dim[1]);
    var intermed_t = new Float64Array(2*dim[0]*dim[1]);
    var fft_trans = new Float64Array(2*dim[0]*dim[1]);

    var fft_type = 'complex';
    var row_len = 2*dim[0];

    for(var i = 0; i < dim[1]; i++){

        fft0.simple(frow,kspace.slice(i*row_len,(i+1)*row_len),fft_type);
        copy_into(frow,intermed,2*dim[0],i*2*dim[0]);
        
    }

    //transpose the row transformed
    transpose_cplx(intermed,dim,intermed_t);

    //fft each of the rows (which are columns from the original image)
    var col_len = 2*dim[1];
    var fcol = new Float64Array(col_len);
    for(var i = 0; i < dim[0]; i++){
        fft1.simple(fcol,intermed_t.slice(i*col_len,(i+1)*col_len),'complex');
        copy_into(fcol,fft_trans,col_len,i*col_len);
    }
    //transpose back to original dimentions
    dim_t = new Float64Array([dim[1], dim[0]]);
    transpose_cplx(fft_trans,dim_t,img_space);
    
    //perform scaling
    scale = 1/(dim[0]*dim[1]);
    for(var i = 0; i < img_space.length; i++){
        img_space[i] = img_space[i]*scale; 
    }
}


//object that stores offset and stretch and transforms coordinate spaces
//offset is applied before stretch
img_proc.CoordTransform = function(oset,stretch) {
    if(oset[0].length != undefined) {
        this.offset = [oset[0][0],oset[1][0],oset[2][0]];
    } else {
        this.offset = oset;
    }

    this.mat = [[stretch[0][0],stretch[0][1],stretch[0][2],this.offset[0]],
                [stretch[1][0],stretch[1][1],stretch[1][2],this.offset[1]],
                [stretch[2][0],stretch[2][1],stretch[2][2],this.offset[2]],
                [0,0,0,1]];

    this.trans = function(xyz){
        var v; 
        if(xyz[0].length != undefined){
            v = [xyz[0][0],xyz[1][0],xyz[2][0],1];
        } else {
            v = [xyz[0],xyz[1],xyz[2],1];
        }
        return numeric.dot(this.mat,v).slice(0,3);
    }

    this.inv = function(xyz){
        var v; 
        if(xyz[0].length != undefined){
            v = [xyz[0][0],xyz[1][0],xyz[2][0],1];
        } else {
            v = [xyz[0],xyz[1],xyz[2],1];
        }
        return numeric.dot(numeric.inv(this.mat),v).slice(0,3);
    }
}


//returns the index of a pixel at coordinate coord when image stack is stored in a 1D vector
img_proc.IndexTransform = function (dim){
    //return (a+dim[0]*b+c*dim[0]*dim[1]);

    this.scale = new Int16Array(dim.length);
    this.scale[0] = 1;
    for(var i = 1; i < dim.length; i++){
        this.scale[i] = this.scale[i-1] * dim[i-1];
    }

    this.get_index = function(coord){
        var idx = 0;
        for(var i = 0; i < coord.length; i++){
            idx += this.scale[i]*coord[i];
        }
        return idx;
    }//this.get_index

    this.get_coord = function(index){
        var coord = new Int16Array(this.scale.length);
        var hold;
        for(var i = coord.length-1; i >= 0; i--){
            hold = index;
            index = hold % this.scale[i];
            coord[i] = (hold-index) / this.scale[i];
        }
        return coord;
    }//this.get_coord
}

//returns values from a continuous n-D gaussian distribution
//using function val()
//function is normalized s.t. val(center) = 1
img_proc.N_Dim_Gauss = function (center,stddevs){
    if(center.length != stddevs.length) { throw 'Error: N_Dim_Gauss parameters must be of same dimension'; }
    this.cent = center;
    this.dev = stddevs;
    this.val = function(coord){
        var sum = 0;
        for(var i = 0; i < this.cent.length; ++i){
            sum += (coord[i]-this.cent[i])*(coord[i]-this.cent[i])/(2*this.dev[i]);
        }
        return Math.exp(-1*sum);
    }
}

img_proc.Plane = function(nrml, p){
    var nrml_len = Math.sqrt(nrml[0]*nrml[0] + nrml[1]*nrml[1] + nrml[2]*nrml[2]);
    this.normal = nrml;
    
    
    
    this.pt = p;
    this.scalar = -1*(nrml[0]*this.pt[0] + nrml[1]*this.pt[1] + nrml[2]*this.pt[2]);

    for(var i = 0; i < 3; i++){
        this.normal[i] /= nrml_len;
    }

    this.scalar /= nrml_len;

    this.distance = function( point ){
        var D = this.normal[0]*point[0] + this.normal[1]*point[1] + this.normal[2]*point[2] + this.scalar;
        return D;
    }

    this.reflect = function ( point ){
        var d  = this.distance(point);
        var refl = new Float32Array(3);
        for(var i = 0 ; i < 3; i++){
            refl[i] = point[i] - (2*this.normal[i]*d);
        }
        return refl;
    }
}

