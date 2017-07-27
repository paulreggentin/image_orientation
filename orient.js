math = require('./kmath.js');
//require("./complex.js");
fs = require("fs");
require('./nifti.js');
numeric = require('./numeric-1.2.6.js')
img_proc = require('./img_processing.js')

var orient = {}

module.exports = orient;


/*
returns information about the orientation of an object in an image
assumes only one object is present
ARGUMENTS:
            img: a float vector representing an image stack
            dim: an integer array representing the dimensions of img
            rw_offset: a 3-element vector representing the real world coordinates of img[0][0][0] 
            rw_dim: a 3x3 matrix that transforms image coordinates into real world coordinates.
            center: a 3-elt vector that will store the center of mass of the image
            normal: a 3-elt vector that will store the normal vector to the symmetry plane of the image
*/
orient.locate = function(img, dim, rw_offset, rw_dim, center, normal){
    var st = performance.now();

    //set up coordinate transforms
    var img_to_rw = new img_proc.CoordTransform(rw_offset,rw_dim);
   


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

    var thresh = min + bin_width * img_proc.otsu(hist,img.length);
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

    var en = performance.now();
    console.log('thresh:',en-st);
    var st = performance.now();

    var tot_mass = 0;
    //var cntr_mass = new Float64Array(3);
    var rwcntr_mass = new Float64Array(3);
//find the center of mass
    
    var idx_cnvrt = new img_proc.IndexTransform(dim);
    for(var c = 0; c < dim[2]; c++){
        for(var b = 0; b < dim[1]; b++){
            for(var a = 0; a < dim[0]; a++){
                var idx = idx_cnvrt.get_index([[a],[b],[c]]);
                if(mask[idx]){
                    tot_mass += img[idx];
                    var r = img_to_rw.trans([a,b,c]);
                    /*
                    cntr_mass[0] += img[idx]*a;
                    cntr_mass[1] += img[idx]*b;
                    cntr_mass[2] += img[idx]*c;
                    */
                    rwcntr_mass[0] += img[idx]*r[0];
                    rwcntr_mass[1] += img[idx]*r[1];
                    rwcntr_mass[2] += img[idx]*r[2];
                }
            }
        }
    }
    /*
    cntr_mass[0] /=tot_mass;
    cntr_mass[1] /=tot_mass;
    cntr_mass[2] /=tot_mass;
    */
    rwcntr_mass[0] /=tot_mass;
    rwcntr_mass[1] /=tot_mass;
    rwcntr_mass[2] /=tot_mass;

    //var rw_cntr_mass2 = img_to_rw.trans(cntr_mass);

    var en = performance.now();
    console.log('center',en-st);
    var st = performance.now();
    
//find the orientation via the principal inertial axes
    var com_offset = new Float64Array([rw_offset[0]-rwcntr_mass[0],rw_offset[1]-rwcntr_mass[1],rw_offset[2]-rwcntr_mass[2]]);
    var img_to_rw_centered = new img_proc.CoordTransform(com_offset,rw_dim);
    I = [[0,0,0],[0,0,0],[0,0,0]];
    for(var c = 0; c < dim[2]; c++){
        for(var b = 0; b < dim[1]; b++){
            for(var a = 0; a < dim[0]; a++){
                var idx = idx_cnvrt.get_index([[a],[b],[c]]);
                if(mask[idx]){
                    var r = img_to_rw_centered.trans([a,b,c]);
                    var aa = r[0];
                    var bb = r[1];
                    var cc = r[2];
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
    console.log('orient',en-st);


    var ex = img.slice(35*slicelen,36*slicelen);
    var scale = ev.E.x[0][axis]/ev.E.x[1][axis];
    for(var x = 0; x < dim[1]; x++){
        var y = math.round((x-rwcntr_mass[1])*scale+rwcntr_mass[0]);
        ex[y+dim[0]*x] = 2*max;
    }
    //img_proc.exportimg(ex,'/home/paul/Documents/inertia_tensor/c.csv');
    
    var princ_axes = [[0,0,0],[0,0,0],[0,0,0]];


    for(var i = 0; i < 3; i++){
        princ_axes[0][i] = ev.E.x[i][0];
        princ_axes[1][i] = ev.E.x[i][1];
        princ_axes[2][i] = ev.E.x[i][2];
    }
    //orientation[0] = rwcntr_mass;
    //orientation[1] = rwcntr_mass[1];
    //orientation[2] = Math.atan(scale);
    var n = img_proc.find_symmetry_plane(img,mask,img_to_rw,idx_cnvrt, princ_axes, center);

    for(var  i= 0; i < 3; i++){
        center[i] = rwcntr_mass[i];
        normal[i] = princ_axes[n][i];
    }
    
}//orient