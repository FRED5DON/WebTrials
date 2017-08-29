var CanvasContext = (function () {
    function CanvasContext(element) {
        this.element = element;
    }
    CanvasContext.prototype.getWebGLContext = function (isDebug) {
        if (isDebug === void 0) { isDebug = false; }
        var canvas = document.querySelector(this.element);
        var supports = ["experimental-webgl", "webgl", "webkit-3d", "moz-webgl"];
        var context;
        for (var i = 0; i < supports.length; i++) {
            context = canvas.getContext(supports[i], isDebug);
            if (context) {
                break;
            }
        }
        context.viewportWidth = canvas.width;
        context.viewportHeight = canvas.height;
        return context;
    };
    //声明初始化着色器函数
    /**
     *
     * @param gl
     * @param vertexShader
     * @param fragmentShader
     * @returns {WebGLProgram}
     */
    CanvasContext.initShader = function (gl, vertexShader, fragmentShader) {
        //把获取到的shader告诉显卡
        //向gl申请一个程序,把shader交给这个程序
        var program = gl.createProgram();
        gl.attachShader(program, vertexShader);
        gl.attachShader(program, fragmentShader);
        //gl链接这个程序
        gl.linkProgram(program);
        //交付给画笔
        gl.useProgram(program);
        return program;
    };
    /**
     *
     * @param gl
     * @param selector
     * @returns {any}
     */
    CanvasContext.getShader = function (gl, selector) {
        var scriptEl = document.querySelector(selector);
        var script = scriptEl.innerText;
        var shader;
        if (scriptEl.getAttribute('type') == 'o-shader/x-vertex') {
            return CanvasContext.getShaderByScripts(gl, script, gl.VERTEX_SHADER);
        }
        else if (scriptEl.getAttribute('type') == 'o-shader/x-fragment') {
            return CanvasContext.getShaderByScripts(gl, script, gl.FRAGMENT_SHADER);
        }
        else {
            return void 0;
        }
    };
    /**
     *
     * @param gl
     * @param script
     * @param type
     * @returns {any}
     */
    CanvasContext.getShaderByScripts = function (gl, script, shaderType) {
        var shader = gl.createShader(shaderType);
        //为gl设置shader、脚本资源
        gl.shaderSource(shader, script);
        gl.compileShader(shader);
        if (gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
            return shader;
        }
        return void 0;
    };
    /**
     *
     * @param gl
     */
    CanvasContext.prototype.draw = function (gl, num) {
        if (num === void 0) { num = 1; }
        gl.viewport(0, 0, gl.viewportWidth, gl.viewportHeight);
        //设置canvas的背景色为黑色
        gl.clearColor(1.0, 1.0, 1.0, 1.0);
        gl.clear(gl.COLOR_BUFFER_BIT);
        gl.drawArrays(gl.POINTS, 0, num);
        //mat4.perspective(null, 45, gl.viewportWidth / gl.viewportHeight, 0.1, 100.0);
    };
    CanvasContext.prototype.drawFragment = function (gl, num) {
        gl.viewport(0, 0, gl.viewportWidth, gl.viewportHeight);
        gl.clear(gl.COLOR_BUFFER_BIT);
        gl.clearDepth(gl.DEPTH_BUFFER_BIT);
        gl.drawArrays(gl.LINE_STRIP, 0, num);
        //mat4.perspective(null, 45, gl.viewportWidth / gl.viewportHeight, 0.1, 100.0);
    };
    return CanvasContext;
})();
//# sourceMappingURL=index.js.map