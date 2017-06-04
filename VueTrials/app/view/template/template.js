;(function (Vue) {

    if (typeof Vue.template != 'object') {
        Vue.template = {};
    }

    Vue.component('app-nav', {
        template: '<nav :class="{\'nav-bar\':nav}"></nav>',
        data: function () {
            return {
                nav: false
            };
        },
        methods: {
            showNav: function (isShow) {
                this.nav = isShow;
            }
        }
    });

    Vue.component('app-content', {
        data: function () {
            return {
                content: ''
            };
        },
        template: '<div class="app-content" >{{getContent}}</div>',
        computed: {
            getContent: function () {
                return this.content.trim();
            }
        },
        methods: {
            setContent: function (outContent) {
                this.content = outContent;
            }
        }
    });

    Vue.component('app-frame', {
        template: '<div :class="stateMask"><slot></slot></div>',
        data: function () {
            return {
                error: false,
                loading: true
            };
        },
        computed: {
            stateMask: function () {
                return {
                    screen: true,
                    loading: this.loading,
                    error: this.error,
                };
            }
        }
    });

    Vue.component('app-tabbar', {
        template: '<div>Bottom</div>'
    });

})(Vue);