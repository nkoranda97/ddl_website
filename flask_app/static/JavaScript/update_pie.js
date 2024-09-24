const data = {...original_data};
const cutoff = slider.value;
const total = data.counts.reduce((a, b) => a + b, 0);

let new_data = {counts: [], angle: [], percentage: [], color: [], [x]: []};
let other_count = 0;

for (let i = 0; i < data.counts.length; i++) {
    const percentage = data.counts[i] / total * 100;
    if (percentage < cutoff) {
        other_count += data.counts[i];
    } else {
        new_data.counts.push(data.counts[i]);
        new_data.angle.push(data.angle[i]);
        new_data.percentage.push(percentage);
        new_data.color.push(data.color[i]);
        new_data[x].push(data[x][i]);
    }
}

if (other_count > 0) {
    new_data.counts.push(other_count);
    new_data.angle.push(other_count / total * 2 * Math.PI);
    new_data.percentage.push(other_count / total * 100);
    new_data.color.push('#d3d3d3');  // Color for "Other" category
    new_data[x].push('Other');
}

source.data = new_data;