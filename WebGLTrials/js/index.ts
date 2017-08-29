class CanvasContext {

    constructor(public element:string) {
    }


    getWebGLContext(isDebug:boolean = false) {
        var canvas = document.querySelector(this.element);
        var supports=["experimental-webgl", "webgl" , "webkit-3d", "moz-webgl"];
        var context;
        for(var i=0;i<supports.length;i++){
            context= canvas.getContext(supports[i], isDebug);
            if(context){
                break;
            }
        }
        context.viewportWidth = canvas.width;
        context.viewportHeight = canvas.height;
        return context;
    }

    //声明初始化着色器函数
    /**
     *
     * @param gl
     * @param vertexShader
     * @param fragmentShader
     * @returns {WebGLProgram}
     */
    public static initShader(gl:WebGLRenderingContext, vertexShader:WebGLShader, fragmentShader:WebGLShader) {
        //把获取到的shader告诉显卡
        //向gl申请一个程序,把shader交给这个程序
        var program:WebGLProgram = gl.createProgram();
        gl.attachShader(program, vertexShader);
        gl.attachShader(program, fragmentShader);
        //gl链接这个程序
        gl.linkProgram(program);
        //交付给画笔
        gl.useProgram(program);
        return program;
    }

    /**
     *
     * @param gl
     * @param selector
     * @returns {any}
     */
    public static getShader(gl:WebGLRenderingContext, selector:string) {
        var scriptEl = document.querySelector(selector);
        var script = scriptEl.innerText;
        var shader:WebGLShader;
        if (scriptEl.getAttribute('type') == 'o-shader/x-vertex') {
            return CanvasContext.getShaderByScripts(gl, script, gl.VERTEX_SHADER);
        } else if (scriptEl.getAttribute('type') == 'o-shader/x-fragment') {
            return CanvasContext.getShaderByScripts(gl, script, gl.FRAGMENT_SHADER);
        } else {
            return void 0;
        }
    }

    /**
     *
     * @param gl
     * @param script
     * @param type
     * @returns {any}
     */
    public static getShaderByScripts(gl:WebGLRenderingContext, script:string, shaderType:number) {
        var shader:WebGLShader = gl.createShader(shaderType);
        //为gl设置shader、脚本资源
        gl.shaderSource(shader, script);
        gl.compileShader(shader);
        if (gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
            return shader;
        }
        return void 0;
    }

    /**
     *
     * @param gl
     */
    draw(gl,num:number=1) {
        gl.viewport(0, 0, gl.viewportWidth, gl.viewportHeight);
        //设置canvas的背景色为黑色
        gl.clearColor(1.0, 1.0, 1.0, 1.0);
        gl.clear(gl.COLOR_BUFFER_BIT);
        gl.drawArrays(gl.POINTS, 0, num);
        //mat4.perspective(null, 45, gl.viewportWidth / gl.viewportHeight, 0.1, 100.0);
    }

    drawFragment(gl,num:number) {
        gl.viewport(0, 0, gl.viewportWidth, gl.viewportHeight);
        gl.clear(gl.COLOR_BUFFER_BIT);
        gl.clearDepth(gl.DEPTH_BUFFER_BIT);
        gl.drawArrays(gl.LINE_STRIP,0, num);

        //mat4.perspective(null, 45, gl.viewportWidth / gl.viewportHeight, 0.1, 100.0);
    }
}

