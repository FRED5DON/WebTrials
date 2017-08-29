/**
 * Created by fred on 2017/8/19.
 */

var M4FHasCSSMatrix = false;
var M4FHasCSSMatrixCopy = false;

if ("WebKitCSSMatrix" in window && ("media" in window && window.media.matchMedium("(-webkit-transform-3d)")) ||
    ("styleMedia" in window && window.styleMedia.matchMedium("(-webkit-transform-3d)"))) {
    M4FHasCSSMatrix = true;
    if ("copy" in WebKitCSSMatrix.prototype)
        M4FHasCSSMatrixCopy = true;
}

class Formulas {

    /**
     * det(Mat2)
     * @param a
     * @param b
     * @param c
     * @param d
     * @returns {number}
     */
    public static determinant2x2(a, b, c, d):number {
        return a * d - b * c;
    }

    /**
     * det(Mat3)
     * @param a1
     * @param a2
     * @param a3
     * @param b1
     * @param b2
     * @param b3
     * @param c1
     * @param c2
     * @param c3
     * @returns {number}
     */
    public static determinant3x3(a1, a2, a3, b1, b2, b3, c1, c2, c3):number {
        return a1 * this.determinant2x2(b2, b3, c2, c3)
            - b1 * this.determinant2x2(a2, a3, c2, c3)
            + c1 * this.determinant2x2(a2, a3, b2, b3);
    }

    /**
     * det(Mat4)
     * @param a1
     * @param a2
     * @param a3
     * @param a4
     * @param b1
     * @param b2
     * @param b3
     * @param b4
     * @param c1
     * @param c2
     * @param c3
     * @param c4
     * @param d1
     * @param d2
     * @param d3
     * @param d4
     * @returns {number}
     */
    public static determinant4x4(a1, a2, a3, a4, b1, b2, b3, b4, c1, c2, c3, c4, d1, d2, d3, d4):number {
        return a1 * this.determinant3x3(b2, b3, b4, c2, c3, c4, d2, d3, d4)
            - b1 * this.determinant3x3(a2, a3, a4, c2, c3, c4, d2, d3, d4)
            + c1 * this.determinant3x3(a2, a3, a4, b2, b3, b4, d2, d3, d4)
            - d1 * this.determinant3x3(a2, a3, a4, b2, b3, b4, c2, c3, c4);
    }


    //array to matrix is has different rule,instead of transpose
    //this arguments is sort by columns
    public static  minors(a1, a2, a3, a4, b1, b2, b3, b4, c1, c2, c3, c4, d1, d2, d3, d4) {
        return new Float32Array([
            this.determinant3x3(b2, b3, b4, c2, c3, c4, d2, d3, d4),
            -this.determinant3x3(b1, b3, b4, c1, c3, c4, d1, d3, d4),
            this.determinant3x3(b1, b2, b4, c1, c2, c4, d1, d2, d4),
            -this.determinant3x3(b1, b2, b3, c1, c2, c3, d1, d2, d3),

            -this.determinant3x3(a2, a3, a4, c2, c3, c4, d2, d3, d4),
            this.determinant3x3(a1, a3, a4, c1, c3, c4, d1, d3, d4),
            -this.determinant3x3(a1, a2, a4, c1, c2, c4, d1, d2, d4),
            this.determinant3x3(a1, a2, a3, c1, c2, c3, d1, d2, d3),

            this.determinant3x3(a2, a3, a4, b2, b3, b4, d2, d3, d4),
            -this.determinant3x3(a1, a3, a4, b1, b3, b4, d1, d3, d4),
            this.determinant3x3(a1, a2, a4, b1, b2, b4, d1, d2, d4),
            -this.determinant3x3(a1, a2, a3, b1, b2, b3, d1, d2, d3),

            -this.determinant3x3(a2, a3, a4, b2, b3, b4, c2, c3, c4),
            this.determinant3x3(a1, a3, a4, b1, b3, b4, c1, c3, c4),
            -this.determinant3x3(a1, a2, a4, b1, b2, b4, c1, c2, c4),
            this.determinant3x3(a1, a2, a3, b1, b2, b3, c1, c2, c3)
        ]);
    }

}

interface Transform {
    scale(x:number, y:number, z:number);

    translate(x:number, y:number, z:number);

    rotate(angle:number, x:number, y:number, z:number);
}

/**
 * Matrix4x4
 * @author fred
 */
interface IMatrix extends Transform {
    /**
     * read a matrix
     * @param mat
     */
    readMatrix(mat:IMatrix):void;

    /**
     * read a array to  matrix
     * @param array
     */
    readArray(array:Float32Array):void;

    /**
     * get a sequence of float array
     * @return Float32Array
     */
    getAsFloat32Array():Float32Array;

    /**
     * make a unit 4x4 matrix　E
     */
    makeIdentity():void;

    /**
     * transfer matrix from row to column or column to row
     * Row Sequence <=> Column Sequence
     */
    transpose():void;

    /**
     * invert the matrix to a revert matrix
     */
    invert():void;

    /**
     * multiply the matrix by the passed frustum values on the right
     * @param left
     * @param right
     * @param bottom
     * @param top
     * @param near
     * @param far
     */
    frustum(left:number, right:number,
            bottom:number, top:number,
            near:number, far:number):void;

    /**
     * multiply the matrix by the passed ortho values on the right
     * @param left
     * @param right
     * @param bottom
     * @param top
     * @param near
     * @param far
     */
    ortho(left:number, right:number,
          bottom:number, top:number,
          near:number, far:number):void;


    lookat(eyex, eyey, eyez, centerx, centery, centerz, upx, upy, upz):void;
    /**
     *
     * @param matL
     * @param matR
     */
    multipleM(matL:IMatrix, matR:IMatrix):IMatrix;

}


class Matrix4 implements IMatrix {

    //raw and operated matrix
    private static SIZE = 4;
    private static ELEMENTS_SIZE = Math.pow(Matrix4.SIZE, 2);
    private $matrix;

    /**
     *
     * @param mat
     */
    constructor(mat:Matrix4 = null) {
        if (M4FHasCSSMatrix) {
            this.$matrix = new WebKitCSSMatrix;
        }
        else {
            this.$matrix = {};
        }
        if (mat) {
            this.readMatrix(mat);
        } else {
            this.makeIdentity();
        }
    }

    readMatrix(mat4:Matrix4):void {
        var matrix = mat4.$matrix;
        if (matrix) {
            this.$matrix.m11 = matrix.m11;
            this.$matrix.m12 = matrix.m12;
            this.$matrix.m13 = matrix.m13;
            this.$matrix.m14 = matrix.m14;

            this.$matrix.m21 = matrix.m21;
            this.$matrix.m22 = matrix.m22;
            this.$matrix.m23 = matrix.m23;
            this.$matrix.m24 = matrix.m24;

            this.$matrix.m31 = matrix.m31;
            this.$matrix.m32 = matrix.m32;
            this.$matrix.m33 = matrix.m33;
            this.$matrix.m34 = matrix.m34;

            this.$matrix.m41 = matrix.m41;
            this.$matrix.m42 = matrix.m42;
            this.$matrix.m43 = matrix.m43;
            this.$matrix.m44 = matrix.m44;
        }
    }

    readArray(array:Float32Array):void {
        if ("length" in array && array.length >= Matrix4.ELEMENTS_SIZE) {
            this.$matrix.m11 = array[0];
            this.$matrix.m12 = array[1];
            this.$matrix.m13 = array[2];
            this.$matrix.m14 = array[3];

            this.$matrix.m21 = array[4];
            this.$matrix.m22 = array[5];
            this.$matrix.m23 = array[6];
            this.$matrix.m24 = array[7];

            this.$matrix.m31 = array[8];
            this.$matrix.m32 = array[9];
            this.$matrix.m33 = array[10];
            this.$matrix.m34 = array[11];

            this.$matrix.m41 = array[12];
            this.$matrix.m42 = array[13];
            this.$matrix.m43 = array[14];
            this.$matrix.m44 = array[15];
        }
    }

    getAsFloat32Array():Float32Array {
        if (M4FHasCSSMatrixCopy) {
            var array = new Float32Array(16);
            this.$matrix.copy(array);
            return array;
        }
        return new Float32Array(this.getAsArray());
    }

    getAsArray() {
        return [
            this.$matrix.m11, this.$matrix.m12, this.$matrix.m13, this.$matrix.m14,
            this.$matrix.m21, this.$matrix.m22, this.$matrix.m23, this.$matrix.m24,
            this.$matrix.m31, this.$matrix.m32, this.$matrix.m33, this.$matrix.m34,
            this.$matrix.m41, this.$matrix.m42, this.$matrix.m43, this.$matrix.m44
        ];
    }

    makeIdentity():void {
        this.$matrix.m11 = 1;
        this.$matrix.m12 = 0;
        this.$matrix.m13 = 0;
        this.$matrix.m14 = 0;

        this.$matrix.m21 = 0;
        this.$matrix.m22 = 1;
        this.$matrix.m23 = 0;
        this.$matrix.m24 = 0;

        this.$matrix.m31 = 0;
        this.$matrix.m32 = 0;
        this.$matrix.m33 = 1;
        this.$matrix.m34 = 0;

        this.$matrix.m41 = 0;
        this.$matrix.m42 = 0;
        this.$matrix.m43 = 0;
        this.$matrix.m44 = 1;
    }

    // row Sequence -> column Sequence
    transpose():void {
        var tmp = this.$matrix.m12;
        this.$matrix.m12 = this.$matrix.m21;
        this.$matrix.m21 = tmp;

        tmp = this.$matrix.m13;
        this.$matrix.m13 = this.$matrix.m31;
        this.$matrix.m31 = tmp;

        tmp = this.$matrix.m14;
        this.$matrix.m14 = this.$matrix.m41;
        this.$matrix.m41 = tmp;

        tmp = this.$matrix.m23;
        this.$matrix.m23 = this.$matrix.m32;
        this.$matrix.m32 = tmp;

        tmp = this.$matrix.m24;
        this.$matrix.m24 = this.$matrix.m42;
        this.$matrix.m42 = tmp;

        tmp = this.$matrix.m34;
        this.$matrix.m34 = this.$matrix.m43;
        this.$matrix.m43 = tmp;
    }

    //invert the matrix
    invert():void {
        if (M4FHasCSSMatrix) {
            this.$matrix = this.$matrix.inverse();
            return;
        }
        //Calculate the 4x4 determinant
        //Cuz only have simple 2x2 or 3x3 determinant formula
        //need convert to calculate sth Multipy with Calculate the 3x3
        let det = this.$determinant();
        //aspect is dot 8 nums
        if (Math.abs(det) < 1e-8)
            return;
        // A x A伴随 = det(A)I
        //即:A逆=A伴随/det(A)
        //(1) 求:A伴随矩阵
        //1.1 Aij为aij的代数余子式 Aij=(-1)的(i+j)次方 x Mij余子式
        //1.2 Mij余子式=在A中除去aij元素所在的行和列组成的n-1阶子式
        this.$minors();
        //(2) /det(A)
        this.$matrix.m11 /= det;
        this.$matrix.m12 /= det;
        this.$matrix.m13 /= det;
        this.$matrix.m14 /= det;

        this.$matrix.m21 /= det;
        this.$matrix.m22 /= det;
        this.$matrix.m23 /= det;
        this.$matrix.m24 /= det;

        this.$matrix.m31 /= det;
        this.$matrix.m32 /= det;
        this.$matrix.m33 /= det;
        this.$matrix.m34 /= det;

        this.$matrix.m41 /= det;
        this.$matrix.m42 /= det;
        this.$matrix.m43 /= det;
        this.$matrix.m44 /= det;
    }


    frustum(left:number, right:number, bottom:number, top:number, near:number, far:number):void {
        var matrix = new Matrix4();
        var A = (right + left) / (right - left);
        var B = (top + bottom) / (top - bottom);
        var C = -(far + near) / (far - near);
        var D = -(2 * far * near) / (far - near);

        matrix.$matrix.m11 = (2 * near) / (right - left);
        matrix.$matrix.m12 = 0;
        matrix.$matrix.m13 = 0;
        matrix.$matrix.m14 = 0;

        matrix.$matrix.m21 = 0;
        matrix.$matrix.m22 = 2 * near / (top - bottom);
        matrix.$matrix.m23 = 0;
        matrix.$matrix.m24 = 0;

        matrix.$matrix.m31 = A;
        matrix.$matrix.m32 = B;
        matrix.$matrix.m33 = C;
        matrix.$matrix.m34 = -1;

        matrix.$matrix.m41 = 0;
        matrix.$matrix.m42 = 0;
        matrix.$matrix.m43 = D;
        matrix.$matrix.m44 = 0;
        this.$multiply(matrix);
    }

    ortho(left:number, right:number, bottom:number, top:number, near:number, far:number) {
        var matrix = new Matrix4();
        matrix.$matrix.m11 = 2 / (right - left);
        matrix.$matrix.m12 = 0;
        matrix.$matrix.m13 = 0;
        matrix.$matrix.m14 = -(right + left) / (right - left);

        matrix.$matrix.m21 = 0;
        matrix.$matrix.m22 = 2 / (top - bottom);
        matrix.$matrix.m23 = 0;
        matrix.$matrix.m24 = -(top + bottom) / (top - bottom);

        matrix.$matrix.m31 = 0;
        matrix.$matrix.m32 = 0;
        matrix.$matrix.m33 = -2 / (far - near);
        matrix.$matrix.m34 = -(far + near) / (far - near);

        matrix.$matrix.m41 = 0;
        matrix.$matrix.m42 = 0;
        matrix.$matrix.m43 = 0;
        matrix.$matrix.m44 = 1;
        this.$multiply(matrix);
    }

    lookat(eyex, eyey, eyez, centerx, centery, centerz, upx, upy, upz) {

        var matrix = new Matrix4();

        // Make rotation matrix
        // Z vector
        var xx, xy, xz;
        var zx = eyex - centerx;
        var zy = eyey - centery;
        var zz = eyez - centerz;
        var mag = Math.sqrt(zx * zx + zy * zy + zz * zz);
        if (mag) {
            zx /= mag;
            zy /= mag;
            zz /= mag;
        }

        // Y vector
        var yx = upx;
        var yy = upy;
        var yz = upz;

        // X vector = Y cross Z
        xx = yy * zz - yz * zy;
        xy = -yx * zz + yz * zx;
        xz = yx * zy - yy * zx;

        // Recompute Y = Z cross X
        yx = zy * xz - zz * xy;
        yy = -zx * xz + zz * xx;
        yz = zx * xy - zy * xx;

        // cross product gives area of parallelogram, which is < 1.0 for
        // non-perpendicular unit-length vectors; so normalize x, y here

        mag = Math.sqrt(xx * xx + xy * xy + xz * xz);
        if (mag) {
            xx /= mag;
            xy /= mag;
            xz /= mag;
        }

        mag = Math.sqrt(yx * yx + yy * yy + yz * yz);
        if (mag) {
            yx /= mag;
            yy /= mag;
            yz /= mag;
        }

        matrix.$matrix.m11 = xx;
        matrix.$matrix.m12 = xy;
        matrix.$matrix.m13 = xz;
        matrix.$matrix.m14 = 0;

        matrix.$matrix.m21 = yx;
        matrix.$matrix.m22 = yy;
        matrix.$matrix.m23 = yz;
        matrix.$matrix.m24 = 0;

        matrix.$matrix.m31 = zx;
        matrix.$matrix.m32 = zy;
        matrix.$matrix.m33 = zz;
        matrix.$matrix.m34 = 0;

        matrix.$matrix.m41 = 0;
        matrix.$matrix.m42 = 0;
        matrix.$matrix.m43 = 0;
        matrix.$matrix.m44 = 1;
        matrix.translate(-eyex, -eyey, -eyez);

        this.$multiply(matrix);
    }

    perspective(fovy, aspect, near, far) {
        var matrix = new Matrix4();
        var tanhf = Math.tan(fovy * 0.5 * Math.PI / 180);
        matrix.$matrix.m11 = 1 / tanhf / aspect;
        matrix.$matrix.m12 = 0;
        matrix.$matrix.m13 = 0;
        matrix.$matrix.m14 = 0;

        matrix.$matrix.m21 = 0;
        matrix.$matrix.m22 = 1 / tanhf;
        matrix.$matrix.m23 = 0;
        matrix.$matrix.m24 = 0;

        matrix.$matrix.m31 = 0;
        matrix.$matrix.m32 = 0;
        matrix.$matrix.m33 = -(far + near ) / (far - near);
        matrix.$matrix.m34 = -1;

        matrix.$matrix.m41 = 0;
        matrix.$matrix.m42 = 0;
        matrix.$matrix.m43 = -(2 * far * near ) / (far - near);
        matrix.$matrix.m44 = 0;
        //this.frustum(left, right, bottom, top, zNear, zFar);
        this.$multiply(matrix);
    }

    multipleM(matL:Matrix4, matR:Matrix4):Matrix4 {
        matR.$multiply(matL);
        return matR;
    }

    scale(x:number, y:number, z:number) {
        if (M4FHasCSSMatrix) {
            this.$matrix = this.$matrix.scale(x, y, z);
            return;
        }

        var matrix = new Matrix4();
        matrix.$matrix.m11 = x;
        matrix.$matrix.m22 = y;
        matrix.$matrix.m33 = z;
        this.$multiply(matrix);
    }

    translate(x:number, y:number, z:number) {
        if (M4FHasCSSMatrix) {
            this.$matrix = this.$matrix.translate(x, y, z);
            return;
        }

        var matrix = new Matrix4();
        matrix.$matrix.m41 = x;
        matrix.$matrix.m42 = y;
        matrix.$matrix.m43 = z;
        this.$multiply(matrix);
    }

    rotate(angle:number, x:number, y:number, z:number) {
        if (M4FHasCSSMatrix) {
            this.$matrix = this.$matrix.rotateAxisAngle(x, y, z, angle);
            return;
        }
        // angles are in degrees. Switch to radians
        var rad = angle / 180 * Math.PI;
        rad /= 2;
        var sinA = Math.sin(rad);
        var cosA = Math.cos(rad);
        var sinA2 = sinA * sinA;

        // normalize
        var len = Math.sqrt(x * x + y * y + z * z);
        if (len == 0) {
            // bad vector, just use something reasonable
            x = 0;
            y = 0;
            z = 1;
        } else if (len != 1) {
            x /= len;
            y /= len;
            z /= len;
        }

        var mat = new Matrix4();

        // optimize case where axis is along major axis
        if (x == 1 && y == 0 && z == 0) {
            mat.$matrix.m11 = 1;
            mat.$matrix.m12 = 0;
            mat.$matrix.m13 = 0;
            mat.$matrix.m21 = 0;
            mat.$matrix.m22 = 1 - 2 * sinA2;
            mat.$matrix.m23 = 2 * sinA * cosA;
            mat.$matrix.m31 = 0;
            mat.$matrix.m32 = -2 * sinA * cosA;
            mat.$matrix.m33 = 1 - 2 * sinA2;
            mat.$matrix.m14 = mat.$matrix.m24 = mat.$matrix.m34 = 0;
            mat.$matrix.m41 = mat.$matrix.m42 = mat.$matrix.m43 = 0;
            mat.$matrix.m44 = 1;
        } else if (x == 0 && y == 1 && z == 0) {
            mat.$matrix.m11 = 1 - 2 * sinA2;
            mat.$matrix.m12 = 0;
            mat.$matrix.m13 = -2 * sinA * cosA;
            mat.$matrix.m21 = 0;
            mat.$matrix.m22 = 1;
            mat.$matrix.m23 = 0;
            mat.$matrix.m31 = 2 * sinA * cosA;
            mat.$matrix.m32 = 0;
            mat.$matrix.m33 = 1 - 2 * sinA2;
            mat.$matrix.m14 = mat.$matrix.m24 = mat.$matrix.m34 = 0;
            mat.$matrix.m41 = mat.$matrix.m42 = mat.$matrix.m43 = 0;
            mat.$matrix.m44 = 1;
        } else if (x == 0 && y == 0 && z == 1) {
            //[
            //    cosB, sinB, 0.0, 0.0,
            //    -sinB, cosB, 0.0, 0.0,
            //    0.0, 0.0, 1.0, 0.0,
            //    0.0, 0.0, 0.0, w
            //]
            mat.$matrix.m11 = 1 - 2 * sinA2;
            mat.$matrix.m12 = 2 * sinA * cosA;
            mat.$matrix.m13 = 0;
            mat.$matrix.m21 = -2 * sinA * cosA;
            mat.$matrix.m22 = 1 - 2 * sinA2;
            mat.$matrix.m23 = 0;
            mat.$matrix.m31 = 0;
            mat.$matrix.m32 = 0;
            mat.$matrix.m33 = 1;
            mat.$matrix.m14 = mat.$matrix.m24 = mat.$matrix.m34 = 0;
            mat.$matrix.m41 = mat.$matrix.m42 = mat.$matrix.m43 = 0;
            mat.$matrix.m44 = 1;
        } else {
            var x2 = x * x;
            var y2 = y * y;
            var z2 = z * z;

            mat.$matrix.m11 = 1 - 2 * (y2 + z2) * sinA2;
            mat.$matrix.m12 = 2 * (x * y * sinA2 + z * sinA * cosA);
            mat.$matrix.m13 = 2 * (x * z * sinA2 - y * sinA * cosA);
            mat.$matrix.m21 = 2 * (y * x * sinA2 - z * sinA * cosA);
            mat.$matrix.m22 = 1 - 2 * (z2 + x2) * sinA2;
            mat.$matrix.m23 = 2 * (y * z * sinA2 + x * sinA * cosA);
            mat.$matrix.m31 = 2 * (z * x * sinA2 + y * sinA * cosA);
            mat.$matrix.m32 = 2 * (z * y * sinA2 - x * sinA * cosA);
            mat.$matrix.m33 = 1 - 2 * (x2 + y2) * sinA2;
            mat.$matrix.m14 = mat.$matrix.m24 = mat.$matrix.m34 = 0;
            mat.$matrix.m41 = mat.$matrix.m42 = mat.$matrix.m43 = 0;
            mat.$matrix.m44 = 1;
        }
        this.$multiply(mat);
    }


    private $determinant():number {
        var a1 = this.$matrix.m11;
        var b1 = this.$matrix.m12;
        var c1 = this.$matrix.m13;
        var d1 = this.$matrix.m14;

        var a2 = this.$matrix.m21;
        var b2 = this.$matrix.m22;
        var c2 = this.$matrix.m23;
        var d2 = this.$matrix.m24;

        var a3 = this.$matrix.m31;
        var b3 = this.$matrix.m32;
        var c3 = this.$matrix.m33;
        var d3 = this.$matrix.m34;

        var a4 = this.$matrix.m41;
        var b4 = this.$matrix.m42;
        var c4 = this.$matrix.m43;
        var d4 = this.$matrix.m44;
        return Formulas.determinant4x4(
            a1, a2, a3, a4,
            b1, b2, b3, b4,
            c1, c2, c3, c4,
            d1, d2, d3, d4);
    }


    private $minors():void {
        var a1 = this.$matrix.m11;
        var b1 = this.$matrix.m12;
        var c1 = this.$matrix.m13;
        var d1 = this.$matrix.m14;

        var a2 = this.$matrix.m21;
        var b2 = this.$matrix.m22;
        var c2 = this.$matrix.m23;
        var d2 = this.$matrix.m24;

        var a3 = this.$matrix.m31;
        var b3 = this.$matrix.m32;
        var c3 = this.$matrix.m33;
        var d3 = this.$matrix.m34;

        var a4 = this.$matrix.m41;
        var b4 = this.$matrix.m42;
        var c4 = this.$matrix.m43;
        var d4 = this.$matrix.m44;
        var arrays = Formulas.minors(a1, a2, a3, a4, b1, b2, b3, b4, c1, c2, c3, c4, d1, d2, d3, d4);
        this.readArray(arrays);
    }

    private $multiply(mat:Matrix4):void {
        if (M4FHasCSSMatrix) {
            this.$matrix = this.$matrix.multiply(mat.$matrix);
            return;
        }

        var m11 = (mat.$matrix.m11 * this.$matrix.m11 + mat.$matrix.m12 * this.$matrix.m21
        + mat.$matrix.m13 * this.$matrix.m31 + mat.$matrix.m14 * this.$matrix.m41);
        var m12 = (mat.$matrix.m11 * this.$matrix.m12 + mat.$matrix.m12 * this.$matrix.m22
        + mat.$matrix.m13 * this.$matrix.m32 + mat.$matrix.m14 * this.$matrix.m42);
        var m13 = (mat.$matrix.m11 * this.$matrix.m13 + mat.$matrix.m12 * this.$matrix.m23
        + mat.$matrix.m13 * this.$matrix.m33 + mat.$matrix.m14 * this.$matrix.m43);
        var m14 = (mat.$matrix.m11 * this.$matrix.m14 + mat.$matrix.m12 * this.$matrix.m24
        + mat.$matrix.m13 * this.$matrix.m34 + mat.$matrix.m14 * this.$matrix.m44);

        var m21 = (mat.$matrix.m21 * this.$matrix.m11 + mat.$matrix.m22 * this.$matrix.m21
        + mat.$matrix.m23 * this.$matrix.m31 + mat.$matrix.m24 * this.$matrix.m41);
        var m22 = (mat.$matrix.m21 * this.$matrix.m12 + mat.$matrix.m22 * this.$matrix.m22
        + mat.$matrix.m23 * this.$matrix.m32 + mat.$matrix.m24 * this.$matrix.m42);
        var m23 = (mat.$matrix.m21 * this.$matrix.m13 + mat.$matrix.m22 * this.$matrix.m23
        + mat.$matrix.m23 * this.$matrix.m33 + mat.$matrix.m24 * this.$matrix.m43);
        var m24 = (mat.$matrix.m21 * this.$matrix.m14 + mat.$matrix.m22 * this.$matrix.m24
        + mat.$matrix.m23 * this.$matrix.m34 + mat.$matrix.m24 * this.$matrix.m44);

        var m31 = (mat.$matrix.m31 * this.$matrix.m11 + mat.$matrix.m32 * this.$matrix.m21
        + mat.$matrix.m33 * this.$matrix.m31 + mat.$matrix.m34 * this.$matrix.m41);
        var m32 = (mat.$matrix.m31 * this.$matrix.m12 + mat.$matrix.m32 * this.$matrix.m22
        + mat.$matrix.m33 * this.$matrix.m32 + mat.$matrix.m34 * this.$matrix.m42);
        var m33 = (mat.$matrix.m31 * this.$matrix.m13 + mat.$matrix.m32 * this.$matrix.m23
        + mat.$matrix.m33 * this.$matrix.m33 + mat.$matrix.m34 * this.$matrix.m43);
        var m34 = (mat.$matrix.m31 * this.$matrix.m14 + mat.$matrix.m32 * this.$matrix.m24
        + mat.$matrix.m33 * this.$matrix.m34 + mat.$matrix.m34 * this.$matrix.m44);

        var m41 = (mat.$matrix.m41 * this.$matrix.m11 + mat.$matrix.m42 * this.$matrix.m21
        + mat.$matrix.m43 * this.$matrix.m31 + mat.$matrix.m44 * this.$matrix.m41);
        var m42 = (mat.$matrix.m41 * this.$matrix.m12 + mat.$matrix.m42 * this.$matrix.m22
        + mat.$matrix.m43 * this.$matrix.m32 + mat.$matrix.m44 * this.$matrix.m42);
        var m43 = (mat.$matrix.m41 * this.$matrix.m13 + mat.$matrix.m42 * this.$matrix.m23
        + mat.$matrix.m43 * this.$matrix.m33 + mat.$matrix.m44 * this.$matrix.m43);
        var m44 = (mat.$matrix.m41 * this.$matrix.m14 + mat.$matrix.m42 * this.$matrix.m24
        + mat.$matrix.m43 * this.$matrix.m34 + mat.$matrix.m44 * this.$matrix.m44);

        this.$matrix.m11 = m11;
        this.$matrix.m12 = m12;
        this.$matrix.m13 = m13;
        this.$matrix.m14 = m14;

        this.$matrix.m21 = m21;
        this.$matrix.m22 = m22;
        this.$matrix.m23 = m23;
        this.$matrix.m24 = m24;

        this.$matrix.m31 = m31;
        this.$matrix.m32 = m32;
        this.$matrix.m33 = m33;
        this.$matrix.m34 = m34;

        this.$matrix.m41 = m41;
        this.$matrix.m42 = m42;
        this.$matrix.m43 = m43;
        this.$matrix.m44 = m44;
    }

    static fillMatrix4x4WithZero(a:Float32Array) {
        const count = 16;
        if (!a) {
            a = new Float32Array(count);
        }
        var i = 0;
        while (i++ < count) {
            if (a[i] === void 0) {
                a[i] = 0;
            }
        }
        return a;
    }


}