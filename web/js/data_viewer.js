// data viewer
// contains common routines for displaying cohomology data onm the website


var data_sources = [];
var data_views = [];
var data_descriptors = [];
var cur_view = 0;
var nparams;

function tparam(l, A, b) {
    let res = b.slice();
    for (let i=0;i<nparams;i++) {
        for (let j=0;j<nparams;j++) {
            res[i] += l[j] * A[i][j];
        }
    }
    return res;
}

function ndrange(ranges) {
    if (ranges.length == 0){
        return [ [] ];
    } else {
        let ret0 = ndrange(ranges.slice(1));
        let ret = [];
        for(let j=ranges[0][0];j<=ranges[0][1];j++) {
            for (const rr of ret0) {
                ret.push([j].concat(rr));
            }
        }
        return ret;
    }
}

function viewselect_change() {
    cur_view = parseInt(viewselect.value);
    let curv = data_views[cur_view];
    viewdetail.innerHTML = curv.description;

    // create new value dictionary
    let dtransp =  Array(nparams);
    for (let i=0;i<nparams;i++){
        dtransp[i]=[];
    }
    let ds = data_sources[0];
    let dd = {};
    for (const l of ds.data) {
        let ll = tparam(l, curv.A, curv.b);
        // console.log(ll)
        for (let i=0;i<nparams;i++) {
            dtransp[i].push(ll[i]);
        }
        dd[ ll ] = l[nparams];
    } 

    console.log(dd);
    // find bounds 
    // let minps = dtransp.map( l  => Math.min(...l));
    // let maxps = dtransp.map( l  => Math.max(...l));
    let psrange = dtransp.map( l  => [Math.min(...l), Math.max(...l)]);
    // console.log(dtransp);
    console.log(psrange);
    let s = "";

    // the last two params are iterated in every table
    let tableps = ndrange(psrange.slice(0,nparams-2));
    let intable = psrange.slice(nparams-2, nparams+1);
    console.log(tableps)
    for (const ps of tableps) {
        
        s += String(ps);
        s += "<table class='datatable'><tr><th></th>";
        for (let i=intable[0][0];i<=intable[0][1];i++){
            s += `<th>${i}</th>`;
        }
        s += "</tr>";
        for (let j=intable[1][0];j<=intable[1][1];j++){
            s += `<tr><td>${j}</td>`;
            for (let i=intable[0][0];i<=intable[0][1];i++){
                // console.log(ps.concat( [i,j] ));
                // console.log(dd[ ps.concat( [i,j] ) ]);
                entry = dd[ ps.concat( [i,j] ) ] || "-";
                s += `<td>${entry}</td>`;
            }
            s += "</tr>";
        }
        s += "</table>";
    } 
    tableview.innerHTML = s;


}

window.onload = (e) => {
    nparams = data_sources[0].data[0].length - 1 // we assume there is at least one data source
    var viewdiv = document.getElementById("dataview");

    // create ui
    let s = "<label for ='viewselect'>Choose view:</label> <select name='viewselect' id='viewselect' onchange='viewselect_change()'>"
    for (let j=0;j< data_views.length;j++)
    {
        s += `<option value=${j}>${data_views[j].name}</option>`;
    }
    s+="</select><div id='viewdetail'></div><div id='tableview'></div>";
    viewdiv.innerHTML = s;

    var viewselect = document.getElementById("viewselect");
    var viewdetail = document.getElementById("viewdetail");
    var tableview = document.getElementById("tableview");

    // fill data for the first time
    viewselect_change();
    
};




