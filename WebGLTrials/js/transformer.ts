class Transformer {

    constructor() {
    }

    /**
     *
     * @param x
     * @param y
     * @param z
     * @returns {Float32Array}
     */
    getTranslateArray(x:number = 0.0, y:number = 0.0, z:number = 0) {
        return new Float32Array([
            1.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0,
            x, y, z, 1.0
        ]);
    }

    /**
     *
     * @param angle
     * @param w
     * @returns {Float32Array}
     */
    getRotateZArray(angle:number = 0, w:number = 1.0) {
        var rad = Math.PI * angle / 180.0;
        var cosB = Math.cos(rad);
        var sinB = Math.sin(rad);
        return new Float32Array([
            cosB, sinB, 0.0, 0.0,
            -sinB, cosB, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, w
        ]);
    }

    /**
     *
     * @param a
     * @param b
     * @param c
     * @returns {Float32Array}
     */
    getScaleArray(a:number = 0.0, b:number = 0.0, c:number = 0) {
        return new Float32Array([
            a, 0.0, 0.0, 0.0,
            0.0, b, 0.0, 0.0,
            0.0, 0.0, c, 0.0,
            0.0, 0.0, 0.0, 1.0
        ]);
    }


    /**
     *
     * @param gl
     * @param u_matrix_pos
     * @param x
     * @param y
     * @param z
     * @param w
     */
    translate(gl:WebGLRenderingContext, u_matrix_pos, x:number = 0.0, y:number = 0.0, z:number = 0) {
        gl.uniformMatrix4fv(u_matrix_pos, false, this.getTranslateArray(x, y, z));
    }


    /**
     *
     * @param gl
     * @param u_matrix_pos
     * @param angle
     */
    rotateZ(gl:WebGLRenderingContext, u_matrix_pos, angle:number = 0.0) {
        gl.uniformMatrix4fv(u_matrix_pos, false, this.getRotateZArray(angle));
    }

    /**
     *
     * @param gl
     * @param u_matrix_pos
     * @param x
     * @param y
     * @param z
     */
    scale(gl:WebGLRenderingContext, u_matrix_pos, x:number = 0.0, y:number = 0.0, z:number = 0) {
        gl.uniformMatrix4fv(u_matrix_pos, false, this.getScaleArray(x, y, z));
    }


    /**
     * notice the sequence
     * left x right
     * @param l
     * @param r
     * @returns {Float32Array}
     */
    matrix4x4Mutipy(l:Float32Array, r:Float32Array) {
        var a = Transformer.fillMatrix4x4WithZero(l);
        var b = Transformer.fillMatrix4x4WithZero(r);
        //Matrix Mutiple
        var res = [
            a[0] * b[0] + a[4] * b[1] + a[8] * b[2] + a[12] * b[3],
            a[0] * b[4] + a[4] * b[5] + a[8] * b[6] + a[12] * b[7],
            a[0] * b[8] + a[4] * b[9] + a[8] * b[10] + a[12] * b[11],
            a[0] * b[12] + a[4] * b[13] + a[8] * b[14] + a[12] * b[15],

            a[1] * b[0] + a[5] * b[1] + a[9] * b[2] + a[13] * b[3],
            a[1] * b[4] + a[5] * b[5] + a[9] * b[6] + a[13] * b[7],
            a[1] * b[8] + a[5] * b[9] + a[9] * b[10] + a[13] * b[11],
            a[1] * b[12] + a[5] * b[13] + a[9] * b[14] + a[13] * b[15],

            a[2] * b[0] + a[6] * b[1] + a[10] * b[2] + a[14] * b[3],
            a[2] * b[4] + a[6] * b[5] + a[10] * b[6] + a[14] * b[7],
            a[2] * b[8] + a[6] * b[9] + a[10] * b[10] + a[14] * b[11],
            a[2] * b[12] + a[6] * b[13] + a[10] * b[14] + a[14] * b[15],

            a[3] * b[0] + a[7] * b[1] + a[11] * b[2] + a[15] * b[3],
            a[3] * b[4] + a[7] * b[5] + a[11] * b[6] + a[15] * b[7],
            a[3] * b[8] + a[7] * b[9] + a[11] * b[10] + a[15] * b[11],
            a[3] * b[12] + a[7] * b[13] + a[11] * b[14] + a[15] * b[15]
        ];
        return new Float32Array([
            res[0], res[4], res[8], res[12],
            res[1], res[5], res[9], res[13],
            res[2], res[6], res[10], res[14],
            res[3], res[7], res[11], res[15]
        ]);
        //return new Float32Array([
        //    a[0] * b[0] + a[1] * b[4] + a[2] * b[8] + a[3] * b[12],
        //    a[0] * b[1] + a[1] * b[5] + a[2] * b[9] + a[3] * b[13],
        //    a[0] * b[2] + a[1] * b[6] + a[2] * b[10] + a[3] * b[14],
        //    a[0] * b[3] + a[1] * b[7] + a[2] * b[11] + a[3] * b[15],
        //
        //    a[4] * b[0] + a[5] * b[4] + a[6] * b[8] + a[7] * b[12],
        //    a[4] * b[1] + a[5] * b[5] + a[6] * b[9] + a[7] * b[13],
        //    a[4] * b[2] + a[5] * b[6] + a[6] * b[10] + a[7] * b[14],
        //    a[4] * b[3] + a[5] * b[7] + a[6] * b[11] + a[7] * b[15],
        //
        //    a[8] * b[0] + a[9] * b[4] + a[10] * b[8] + a[11] * b[12],
        //    a[8] * b[1] + a[9] * b[5] + a[10] * b[9] + a[11] * b[13],
        //    a[8] * b[2] + a[9] * b[6] + a[10] * b[10] + a[11] * b[14],
        //    a[8] * b[3] + a[9] * b[7] + a[10] * b[11] + a[11] * b[15],
        //
        //    a[12] * b[0] + a[13] * b[4] + a[14] * b[8] + a[15] * b[12],
        //    a[12] * b[1] + a[13] * b[5] + a[14] * b[9] + a[15] * b[13],
        //    a[12] * b[2] + a[13] * b[6] + a[14] * b[10] + a[15] * b[14],
        //    a[12] * b[3] + a[13] * b[7] + a[14] * b[11] + a[15] * b[15]
        //]);
    }

    /**
     *
     * @param gl
     * @param u_matrix_pos
     * @param matrix
     */
    transform(gl:WebGLRenderingContext, u_matrix_pos, matrix:Float32Array) {
        gl.uniformMatrix4fv(u_matrix_pos, false, matrix);
    }


    /**
     *
     * @param a
     * @returns {Float32Array}
     */
    static fillMatrix4x4WithZero(a:Float32Array) {
        const count = 16;
        var i = 0;
        while (i++ < count) {
            if (a[i] === void 0) {
                a[i] = 0;
            }
        }
        return a;
    }


    perspective(fovy, aspect, zNear, zFar) {
        var top = Math.tan(fovy * Math.PI / 360) * zNear;
        var bottom = -top;
        var left = aspect * bottom;
        var right = aspect * top;
        this.frustum(left, right, bottom, top, zNear, zFar);
    }

    frustum (left, right, bottom, top, near, far)
    {
        var matrix = new J3DIMatrix4();
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

        this.multiply(matrix);
    }

}