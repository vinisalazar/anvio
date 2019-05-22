/**
 *  Canvas drawing functions.
 *
 *  Author: Özcan Esen <ozcanesen@gmail.com>
 *  Copyright 2019, The anvio Project
 *
 * This file is part of anvi'o (<https://github.com/merenlab/anvio>).
 * 
 * Anvi'o is a free software. You can redistribute this program
 * and/or modify it under the terms of the GNU General Public 
 * License as published by the Free Software Foundation, either 
 * version 3 of the License, or (at your option) any later version.
 * 
 * You should have received a copy of the GNU General Public License
 * along with anvi'o. If not, see <http://opensource.org/licenses/GPL-3.0>.
 *
 * @license GPL-3.0+ <http://opensource.org/licenses/GPL-3.0>
 */


//--------------------------------------------------------------------------------------------------
function createCanvas(parent, canvas_id, x, y, height, width) {
    let svgObject = document.getElementById(parent);

    let foreignObject  = document.createElement('foreignObject');
    foreignObject.setAttribute('x', x);
    foreignObject.setAttribute('y', y);
    foreignObject.setAttribute('width', width);
    foreignObject.setAttribute('height', height);

    let canvasBody = document.createElement('xhtml:body');
    canvasBody.setAttribute('width', width);
    canvasBody.setAttribute('height', height);
    canvasBody.setAttribute('syle', 'margin: 0px; padding: 0px; background-color: red;');
    
    let canvas = document.createElement('canvas');
    canvas.setAttribute('width', width);
    canvas.setAttribute('height', height);

    canvasBody.appendChild(canvas);
    foreignObject.appendChild(canvasBody);
    svgObject.appendChild(foreignObject);
}

