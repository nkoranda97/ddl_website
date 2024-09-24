const selected = cb_obj.item;
const data = df_dict;
const group = {};

for (let i = 0; i < data[selected].length; i++) {
    const key = data[selected][i];
    if (group[key] === undefined) {
        group[key] = 0;
    }
    group[key] += 1;
}

const keys = Object.keys(group);
keys.sort((a, b) => group[b] - group[a]);

const counts = keys.map(key => group[key]);
const factors = keys.map(String);

source.data = {
    [selected]: factors,
    counts: counts
};

p.y_range.factors = factors;
p.renderers[0].glyph.y = { field: selected };

p.title.text = selected;
dropdown.label = selected;  // Update the label of the dropdown

// Update hover tooltips
hover.tooltips = [
    ["Category", `@{${selected}}`],
    ["Count", "@counts"]
];