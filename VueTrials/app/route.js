;
const BASE_PATH = 'res/activities';
const SUFFIX = '.html';
const ROUTE_CONFIG = {
    home: '/home'


};
(function (loc) {
    if (typeof urlTools != 'undefined') {
        var p = urlTools.getQueryString('p');
        if (ROUTE_CONFIG[p]) {
            loc.replace(BASE_PATH + ROUTE_CONFIG[p] + SUFFIX);
        }
    }
})(location);