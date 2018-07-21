


var getUrlParameter = function getUrlParameter(sParam) {
    var sPageURL = decodeURIComponent(window.location.search.substring(1)),
        sURLVariables = sPageURL.split('&'),
        sParameterName,
        i;

    for (i = 0; i < sURLVariables.length; i++) {
        sParameterName = sURLVariables[i].split('=');

        if (sParameterName[0] === sParam) {
            return sParameterName[1] === undefined ? true : sParameterName[1];
        }
    }
};

var dedata;

var plotLayout = {
    yaxis: {
	type: 'category',
	//showgrid: true,
	//title: 'samples',
    },
    xaxis: {
	type: 'category',
	//showgrid: true,
	//title: 'cell types',
    },
    margin: {
	l: 150,
	r: 50,
	b: 100,
	t: 50
    },
    showlegend: false,
    hovermode: 'closest'
};

var plotPanel = Ext.create('Ext.panel.Panel', {
    minHeight: 200,
    id: 'plots',
    minWidth: 200,
    gene: 'GAPDH',
    minMarkerSize: 10,
    height: '100%',
    layout: 'fit',
    dockedItems: [{
	xtype: 'toolbar',
	dock: 'top',
	items: [
	    {
		xtype: 'combobox',
		//fieldLabel: 'Gene',
		id: 'geneField',
		emptyText: 'gene name',
		displayField: 'gene',
		queryMode: 'local',
		typeAhead: true,
		minChars: 1,
		selectOnFocus: true,
		listeners: {
		    select: function(el,n,o) {
			plotPanel.showGene(n[0].data.gene);
		    }
		}
	    },{
		text: 'Show gene',
		listeners: {
		    click: function() {
			var gene=Ext.getCmp("geneField").getValue();
			plotPanel.showGene(gene);
		    }
		}
	    },{
		
	    }
	]
    }],
    //html:['<div id="plots"/>' ],
    render: function() {
	var trace1={ 
	    mode: 'markers',
	    marker: {
		sizemode: 'area'
	    }
	}
	trace1['x']=dedata.cellCounts.type;
	trace1['y']=dedata.cellCounts.sample;
	trace1.marker['size']=dedata.cellCounts.cells;
	// normalize cell to the size of each sample
	var tc={};
	for(i=0;i<trace1['y'].length;i++) {
	    var sn=trace1['y'][i]; 
	    if(!tc.hasOwnProperty(sn)) { 
		tc[sn]=trace1.marker['size'][i]; 
	    } else {
		tc[sn]=tc[sn] + trace1.marker['size'][i]; 
	    }
	}
	var scale=Math.min(plotPanel.getSize().width,plotPanel.getSize().height)/Object.keys(tc).length/2;
	for(i=0;i<trace1['y'].length;i++) {
	    trace1.marker['size'][i]=plotPanel.minMarkerSize+ trace1.marker['size'][i]/tc[ trace1['y'][i] ] * (Object.keys(tc).length)*scale
	}
	// color according to gene counts
	trace1.marker['color']=plotPanel.geneColor(plotPanel.gene);
	plotPanel.plotTrace=trace1;
	plotLayout.yaxis['categoryorder']='array';
	plotLayout.yaxis['categoryarray']=dedata.sampleNames;
	Plotly.newPlot('plots',[trace1],plotLayout)
	plotPanel.on('resize',function(cmp,width,height,oldWidth,oldHeight,opts) {
	    Plotly.Plots.resize('plots');
	})
    },
    showGene: function(gene) {
	plotPanel.gene=gene;
	plotPanel.plotTrace.marker['color']=plotPanel.geneColor(plotPanel.gene);
	Plotly.newPlot('plots',[plotPanel.plotTrace],plotLayout)
    },
    geneColor: function(gene) {
	var val=dedata.geneCounts[gene];
	//for(i=0;i<val.length;i++) { val[i]=Math.log(val[i]+1) }
	// get quantiles
	var trim=0.01;
	var qrng=math.quantileSeq(val,[1.0-trim,1])
	var vmean=math.mean(val);
	var color=d3.scale.linear().domain([0,vmean,qrng[0],qrng[1]])
	    .interpolate(d3.interpolateLab)
	    .range(['blue','#AAAAAA','red','red']);
	var cols=[];
	for(i=0;i<val.length;i++) {
	    cols[i]=color(val[i]);
	}
	return(cols);
    }
    
    
})

Ext.application({
    name: 'Fiddle',
    
    launch: function() {
	

	
	var viewport = Ext.create('Ext.Viewport', {
            layout: {
		type: 'border',
		padding: 0
            },
            defaults: {
		split: true
            },
            items: [{ region: 'center',
		      xtype: 'panel',
		      layout: 'fit',
		      items: [plotPanel]
		    }]
	}).show();
	
	// load data
	plotPanel.setLoading(true);
	Ext.Ajax.request({
	    url: getUrlParameter('d'),
	    headers: { "Accept-Encoding" : "gzip" },
	    success: function(response, opts) {
		dedata = Ext.decode(response.responseText);
		// show the initial plot
		//console.log(dedata.cellCounts);
		Ext.define('GeneNames', {
		    extend: 'Ext.data.Model',
		    fields: [
			{name: 'gene', convert: function(value, record) {
			    return record.raw;
			}}
		    ]
		});
		var gstore=Ext.create('Ext.data.ArrayStore',{
		    autoload: true,
		    fields: ['gene'],
		    model: 'GeneNames',
		    data: Object.keys(dedata.geneCounts)
		})
		Ext.getCmp("geneField").bindStore(gstore);
		plotPanel.setLoading(false);
		plotPanel.render();
	    },
	    failure: function(response, opts) {
		console.log('server-side failure with status code ' + response.status);
		plotPanel.setLoading(false);
	    }
	});
	

    }
});
