/**
 * The Unmatrix jQuery function with the sylvester matrix transformation library included
 * All dead sylvester code has been deleted
 */
;(function ($, window, document) {
    var Sylvester = {
        version: "0.1.3",
        precision: 1e-6
    };

    function Vector() {}
    Vector.prototype = {

        // Returns the modulus ('length') of the vector
        modulus: function() {
            return Math.sqrt(this.dot(this));
        },

        // Returns a copy of the vector
        dup: function() {
            return Vector.create(this.elements);
        },

        // Maps the vector to another vector according to the given function
        map: function(fn) {
            var elements = [];
            this.each(function(x, i) {
                elements.push(fn(x, i));
            });
            return Vector.create(elements);
        },

        // Calls the iterator for each element of the vector in turn
        each: function(fn) {
            var n = this.elements.length, k = n, i;
            do { i = k - n;
                fn(this.elements[i], i+1);
            } while (--n);
        },

        // Returns a new vector created by normalizing the receiver
        toUnitVector: function() {
            var r = this.modulus();
            if (r === 0) { return this.dup(); }
            return this.map(function(x) { return x/r; });
        },

        // Returns the result of multiplying the elements of the vector by the argument
        multiply: function(k) {
            return this.map(function(x) { return x*k; });
        },

        // Returns the scalar product of the vector with the argument
        // Both vectors must have equal dimensionality
        dot: function(vector) {
            var V = vector.elements || vector;
            var i, product = 0, n = this.elements.length;
            if (n != V.length) { return null; }
            do { product += this.elements[n-1] * V[n-1]; } while (--n);
            return product;
        },

        // Returns the vector product of the vector with the argument
        // Both vectors must have dimensionality 3
        cross: function(vector) {
            var B = vector.elements || vector;
            if (this.elements.length != 3 || B.length != 3) { return null; }
            var A = this.elements;
            return Vector.create([
                (A[1] * B[2]) - (A[2] * B[1]),
                (A[2] * B[0]) - (A[0] * B[2]),
                (A[0] * B[1]) - (A[1] * B[0])
            ]);
        },

        // Set vector's elements from an array
        setElements: function(els) {
            this.elements = (els.elements || els).slice();
            return this;
        }
    };

    // Constructor function
    Vector.create = function(elements) {
        var V = new Vector();
        return V.setElements(elements);
    };

    function Matrix() {}
    Matrix.prototype = {

        // Returns column k of the matrix as a vector
        col: function(j) {
            if (j > this.elements[0].length) { return null; }
            var col = [], n = this.elements.length, k = n, i;
            do { i = k - n;
                col.push(this.elements[i][j-1]);
            } while (--n);
            return Vector.create(col);
        },

        // Returns a copy of the matrix
        dup: function() {
            return Matrix.create(this.elements);
        },

        // Maps the matrix to another matrix (of the same dimensions) according to the given function
        map: function(fn) {
            var els = [], ni = this.elements.length, ki = ni, i, nj, kj = this.elements[0].length, j;
            do { i = ki - ni;
                nj = kj;
                els[i] = [];
                do { j = kj - nj;
                    els[i][j] = fn(this.elements[i][j], i + 1, j + 1);
                } while (--nj);
            } while (--ni);
            return Matrix.create(els);
        },

        // Returns true iff the matrix can multiply the argument from the left
        canMultiplyFromLeft: function(matrix) {
            var M = matrix.elements || matrix;
            if (typeof(M[0][0]) == 'undefined') { M = Matrix.create(M).elements; }
            // this.columns should equal matrix.rows
            return (this.elements[0].length == M.length);
        },

        // Returns the result of multiplying the matrix from the right by the argument.
        // If the argument is a scalar then just multiply all the elements. If the argument is
        // a vector, a vector is returned, which saves you having to remember calling
        // col(1) on the result.
        multiply: function(matrix) {
            if (!matrix.elements) {
                return this.map(function(x) { return x * matrix; });
            }
            var returnVector = matrix.modulus ? true : false;
            var M = matrix.elements || matrix;
            if (typeof(M[0][0]) == 'undefined') { M = Matrix.create(M).elements; }
            if (!this.canMultiplyFromLeft(M)) { return null; }
            var ni = this.elements.length, ki = ni, i, nj, kj = M[0].length, j;
            var cols = this.elements[0].length, elements = [], sum, nc, c;
            do { i = ki - ni;
                elements[i] = [];
                nj = kj;
                do { j = kj - nj;
                    sum = 0;
                    nc = cols;
                    do { c = cols - nc;
                        sum += this.elements[i][c] * M[c][j];
                    } while (--nc);
                    elements[i][j] = sum;
                } while (--nj);
            } while (--ni);
            var M = Matrix.create(elements);
            return returnVector ? M.col(1) : M;
        },

        // Returns the transpose of the matrix
        transpose: function() {
            var rows = this.elements.length, cols = this.elements[0].length;
            var elements = [], ni = cols, i, nj, j;
            do { i = cols - ni;
                elements[i] = [];
                nj = rows;
                do { j = rows - nj;
                    elements[i][j] = this.elements[j][i];
                } while (--nj);
            } while (--ni);
            return Matrix.create(elements);
        },

        // Returns true iff the matrix is square
        isSquare: function() {
            return (this.elements.length == this.elements[0].length);
        },

        // Make the matrix upper (right) triangular by Gaussian elimination.
        // This method only adds multiples of rows to other rows. No rows are
        // scaled up or switched, and the determinant is preserved.
        toRightTriangular: function() {
            var M = this.dup(), els;
            var n = this.elements.length, k = n, i, np, kp = this.elements[0].length, p;
            do { i = k - n;
                if (M.elements[i][i] == 0) {
                    for (j = i + 1; j < k; j++) {
                        if (M.elements[j][i] != 0) {
                            els = []; np = kp;
                            do { p = kp - np;
                                els.push(M.elements[i][p] + M.elements[j][p]);
                            } while (--np);
                            M.elements[i] = els;
                            break;
                        }
                    }
                }
                if (M.elements[i][i] != 0) {
                    for (j = i + 1; j < k; j++) {
                        var multiplier = M.elements[j][i] / M.elements[i][i];
                        els = []; np = kp;
                        do { p = kp - np;
                            // Elements with column numbers up to an including the number
                            // of the row that we're subtracting can safely be set straight to
                            // zero, since that's the point of this routine and it avoids having
                            // to loop over and correct rounding errors later
                            els.push(p <= i ? 0 : M.elements[j][p] - M.elements[i][p] * multiplier);
                        } while (--np);
                        M.elements[j] = els;
                    }
                }
            } while (--n);
            return M;
        },

        // Returns the determinant for square matrices
        determinant: function() {
            if (!this.isSquare()) { return null; }
            var M = this.toRightTriangular();
            var det = M.elements[0][0], n = M.elements.length - 1, k = n, i;
            do { i = k - n + 1;
                det = det * M.elements[i][i];
            } while (--n);
            return det;
        },

        // Returns true iff the matrix is singular
        isSingular: function() {
            return (this.isSquare() && this.determinant() === 0);
        },

        // Returns the result of attaching the given argument to the right-hand side of the matrix
        augment: function(matrix) {
            var M = matrix.elements || matrix;
            if (typeof(M[0][0]) == 'undefined') { M = Matrix.create(M).elements; }
            var T = this.dup(), cols = T.elements[0].length;
            var ni = T.elements.length, ki = ni, i, nj, kj = M[0].length, j;
            if (ni != M.length) { return null; }
            do { i = ki - ni;
                nj = kj;
                do { j = kj - nj;
                    T.elements[i][cols + j] = M[i][j];
                } while (--nj);
            } while (--ni);
            return T;
        },

        // Returns the inverse (if one exists) using Gauss-Jordan
        inverse: function() {
            if (!this.isSquare() || this.isSingular()) { return null; }
            var ni = this.elements.length, ki = ni, i, j;
            var M = this.augment(Matrix.I(ni)).toRightTriangular();
            var np, kp = M.elements[0].length, p, els, divisor;
            var inverse_elements = [], new_element;
            // Matrix is non-singular so there will be no zeros on the diagonal
            // Cycle through rows from last to first
            do { i = ni - 1;
                // First, normalise diagonal elements to 1
                els = []; np = kp;
                inverse_elements[i] = [];
                divisor = M.elements[i][i];
                do { p = kp - np;
                    new_element = M.elements[i][p] / divisor;
                    els.push(new_element);
                    // Shuffle of the current row of the right hand side into the results
                    // array as it will not be modified by later runs through this loop
                    if (p >= ki) { inverse_elements[i].push(new_element); }
                } while (--np);
                M.elements[i] = els;
                // Then, subtract this row from those above it to
                // give the identity matrix on the left hand side
                for (j = 0; j < i; j++) {
                    els = []; np = kp;
                    do { p = kp - np;
                        els.push(M.elements[j][p] - M.elements[i][p] * M.elements[j][i]);
                    } while (--np);
                    M.elements[j] = els;
                }
            } while (--ni);
            return Matrix.create(inverse_elements);
        },

        // Set the matrix's elements from an array. If the argument passed
        // is a vector, the resulting matrix will be a single column.
        setElements: function(els) {
            var i, elements = els.elements || els;
            if (typeof(elements[0][0]) != 'undefined') {
                var ni = elements.length, ki = ni, nj, kj, j;
                this.elements = [];
                do { i = ki - ni;
                    nj = elements[i].length; kj = nj;
                    this.elements[i] = [];
                    do { j = kj - nj;
                        this.elements[i][j] = elements[i][j];
                    } while (--nj);
                } while(--ni);
                return this;
            }
            var n = elements.length, k = n;
            this.elements = [];
            do { i = k - n;
                this.elements.push([elements[i]]);
            } while (--n);
            return this;
        }
    };

    // Constructor function
    Matrix.create = function(elements) {
        var M = new Matrix();
        return M.setElements(elements);
    };

    // Identity matrix of size n
    Matrix.I = function(n) {
        var els = [], k = n, i, nj, j;
        do { i = k - n;
            els[i] = []; nj = k;
            do { j = k - nj;
                els[i][j] = (i == j) ? 1 : 0;
            } while (--nj);
        } while (--n);
        return Matrix.create(els);
    };

    // Add methods to Sylvester object.
    var MatrixInstance = Matrix.create;
    var VectorInstance = Vector.create;

    // Converts radians to degrees
    var deg = function (rad) {
        return rad * (180 / Math.PI);
    };

    // Returns the determinant of matrix 
    var determinant = function (matrix) {
        return MatrixInstance(matrix).determinant();
    };

    // Returns the inverse of matrix
    var inverse = function (matrix) {
        return MatrixInstance(matrix).inverse().elements;
    };      

    // Returns the transpose of matrix
    var transpose = function (matrix) {
        return MatrixInstance(matrix).transpose().elements;
    };

    // Multiplies vector by matrix and returns the transformed vector
    var multiplyVectorMatrix = function (vector, matrix) {
        return MatrixInstance(matrix).multiply(vector).elements;
    };

    // Returns the length of vector
    var length = function (vector) {
        return VectorInstance(vector).modulus();
    };

    // Normalizes the length of vector to 1
    var normalize = function (vector) {
        return VectorInstance(vector).toUnitVector().elements;
    };

    // Returns the dot product of two points
    var dot = function (vector1, vector2) {
        return VectorInstance(vector1).dot(vector2);
    };

    // Returns the cross product of two vectors
    var cross = function (vector1, vector2) {
        return VectorInstance(vector1).cross(vector2).elements;
    };

    var combine = function (a, b, ascl, bscl) {
        var result = [];
        result[0] = (ascl * a[0]) + (bscl * b[0]);
        result[1] = (ascl * a[1]) + (bscl * b[1]);
        // Both vectors are 3d. Return a 3d vector
        if (a.length === 3 && b.length === 3) {
            result[2] = (ascl * a[2]) + (bscl * b[2]);
        }
        return result;
    };

    // Returns null if the matrix cannot be decomposed, an object if it can
    var unmatrix = function (matrix) {
        var rotateX;
        var rotateY;
        var rotateZ; 
        var scaleX;
        var scaleY;
        var scaleZ;
        var skew;
        var skewX;
        var skewY;
        var translateX; 
        var translateY;
        var translateZ;

        // Normalize the matrix
        if (matrix[3][3] === 0) {
            return null;
        }

        for (var i = 0; i < 4; i++) {
            for (var j = 0; j < 4; j++) {
                matrix[i][j] /= matrix[3][3];
            }
        }

        // perspectiveMatrix is used to solve for perspective, but it also 
        // provides an easy way to test for singularity of the upper 3x3 
        // component
        var perspectiveMatrix = matrix;

        for (var i = 0; i < 3; i++) {
            perspectiveMatrix[i][3] = 0;
        }

        perspectiveMatrix[3][3] = 1;

        if (determinant(perspectiveMatrix) === 0) {
           return null;
        }

        // First, isolate perspective
        var perspective;
        if (matrix[0][3] !== 0 || matrix[1][3] !== 0 || matrix[2][3] !== 0) {
            // rightHandSide is the right hand side of the equation
            var rightHandSide = [];
            rightHandSide[0] = matrix[0][3];
            rightHandSide[1] = matrix[1][3];
            rightHandSide[2] = matrix[2][3];
            rightHandSide[3] = matrix[3][3];

            // Solve the equation by inverting perspectiveMatrix and multiplying 
            // rightHandSide by the inverse
            var inversePerspectiveMatrix = inverse(perspectiveMatrix);
            var transposedInversePerspectiveMatrix = transpose(inversePerspectiveMatrix);
            perspective = multiplyVectorMatrix(rightHandSide, transposedInversePerspectiveMatrix);

            // Clear the perspective partition
            matrix[0][3] = matrix[1][3] = matrix[2][3] = 0;
            matrix[3][3] = 1;
        } else {
            // No perspective
            perspective = [];
            perspective[0] = perspective[1] = perspective[2] = 0;
            perspective[3] = 1;
        }

        // Next take care of translation
        translateX = matrix[3][0];
        //matrix[3][0] = 0;
        translateY = matrix[3][1];
        //matrix[3][1] = 0;
        translateZ = matrix[3][2];
        //matrix[3][2] = 0;

        // Now get scale and shear. "row" is a 3 element array of 3 component 
        // vectors
        var row = [[], [], []];

        for (var i = 0; i < 3; i++) {
            row[i][0] = matrix[i][0];
            row[i][1] = matrix[i][1];
            row[i][2] = matrix[i][2];
        }

        // Compute X scale factor and normalize first row
        scaleX = length(row[0]);
        row[0] = normalize(row[0]);

        // Compute XY shear factor and make 2nd row orthogonal to 1st
        skew = dot(row[0], row[1]);
        row[1] = combine(row[1], row[0], 1.0, -skew);

        // Now, compute Y scale and normalize 2nd row
        scaleY = length(row[1]);
        row[1] = normalize(row[1]);
        skew /= scaleY;

        // Compute XZ and YZ shears, orthogonalize 3rd row
        skewX = dot(row[0], row[2]);
        row[2] = combine(row[2], row[0], 1.0, -skewX);
        skewY = dot(row[1], row[2]);
        row[2] = combine(row[2], row[1], 1.0, -skewY);

        // Next, get Z scale and normalize 3rd row
        scaleZ = length(row[2]);
        row[2] = normalize(row[2]);
        skewX /= scaleZ;
        skewY /= scaleZ;

        // At this point, the matrix (in rows) is orthonormal. Check for a 
        // coordinate system flip. If the determinant is -1, then negate the 
        // matrix and the scaling factors
        var pdum3 = cross(row[1], row[2]);

        if (dot(row[0], pdum3) < 0) {
            for (var i = 0; i < 3; i++) {
                scaleX *= -1;
                row[i][0] *= -1;
                row[i][1] *= -1;
                row[i][2] *= -1;
            }
        }

        // Get the rotations
        rotateY = Math.asin(-row[0][2]);
        if (Math.cos(rotateY) !== 0) {
            rotateX = Math.atan2(row[1][2], row[2][2]);
            rotateZ = Math.atan2(row[0][1], row[0][0]);
        } else {
            rotateX = Math.atan2(-row[2][0], row[1][1]);
            rotateZ = 0;
        }

        return {
            rotate:     deg(rotateZ),
            rotateX:    deg(rotateX),
            rotateY:    deg(rotateY), 
            rotateZ:    deg(rotateZ), 
            scaleX:     scaleX,
            scaleY:     scaleY, 
            scaleZ:     scaleZ,
            skew:       deg(skew),
            skewX:      deg(skewX), 
            skewY:      deg(skewY),
            translateX: translateX, 
            translateY: translateY,
            translateZ: translateZ
        };
    };

    // Returns an object with transform properties
    var getTransform = function (cssTransform) {
        // Check if the transform is 3d
        var is3d = cssTransform.indexOf("matrix3d") > -1;

        // Convert matrix values to an array
        cssTransform = cssTransform.match(/[\d\s\.,\-]+(?!d)/)[0];
        var values = cssTransform.split(",");

        // Convert values to floats
        for (var i = 0, l = values.length; i < l; i++) {
            values[i] = parseFloat(values[i]).toFixed(2);
        }

        // Matrix columns become arrays
        var matrix = is3d ? [ // Create 4x4 3d matrix
                                [values[0],   values[1],  values[2],  values[3]],
                                [values[4],   values[5],  values[6],  values[7]],
                                [values[8],   values[9],  values[10], values[11]],
                                [values[12],  values[13], values[14], values[15]]
                            ] : 
                            [ // Create 4x4 2d matrix
                                [values[0],   values[1],  0,          0],
                                [values[2],   values[3],  0,          0],
                                [0,           0,          1,          0],
                                [values[4],   values[5],  0,          1]
                            ];

        return unmatrix(matrix);
    };

    $.fn.unmatrix = function () {
        var properties = ["transform", "MozTransform", "msTransform", "OTransform", "webkitTransform"];
        var transformProperty = "";

        // Check if browser supports transforms
        $.each(properties, function (index, property) {
            if (property in document.body.style) {
                transformProperty = property;
                return false;
            }
        });

        if (!transformProperty) {
            // Browser does not support transforms
            return false;
        } else {
            // Browser supports transforms
            var transforms = [];

            this.each(function (index, element) {
                var cssTransform = $(element).css(transformProperty);
                var transform = cssTransform !== "none" ? 
                                getTransform(cssTransform) : 
                                {};
                transforms.push(transform);
            });
            return transforms;
        }
    };
})(jQuery, window, document);