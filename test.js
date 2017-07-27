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
numeric = require('./numeric-1.2.6.js');
orient = require('./orient.js');
img_proc = require('./img_processing.js');


var filename = 'rotated.nii'

var raw = fs.readFileSync(filename);

var nii = parse(raw);



center = new Float32Array(3);
normal = new Float32Array(3);


orient.locate(nii.data, nii.sizes, nii.spaceOrigin, nii.spaceDirections, center, normal);

//orient(nii.data, nii.sizes, [5,10,15], [[2,0,0],[0,4,0],[0,0,8]], orientation);

console.log(center,normal);


/*
var sm = new Float32Array(nii.data.length/64);



img_proc.downsample(nii.data,nii.sizes, [4,4,4], sm);

img_proc.exportimg(sm,'sm.csv');
*/

/*
var p = new img_proc.Plane([2,-3,-4],[0,22,0]);

console.log(p.reflect([1,-2,4]));
*/