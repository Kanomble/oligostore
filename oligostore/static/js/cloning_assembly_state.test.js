const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const test = require("node:test");

const selection = require("./cloning_assembly_state");

test("createState stores restriction-ligation form selections", () => {
  const state = selection.createState({
    vector_asset: "sequence_file:1:vec1",
    insert_asset: "sequence_file:2:ins1",
    assembly_strategy: "restriction_ligation",
    is_circular: "1",
    left_enzyme: " EcoRI ",
    right_enzyme: "BamHI",
    selectedEnzymes: ["EcoRI", "EcoRI", "", "BamHI"],
    vector_fragment_index: "2",
    insert_fragment_index: "1",
  });

  assert.equal(state.vectorAsset, "sequence_file:1:vec1");
  assert.equal(state.insertAsset, "sequence_file:2:ins1");
  assert.equal(state.assemblyStrategy, "restriction_ligation");
  assert.equal(state.isCircular, "1");
  assert.equal(state.leftEnzyme, "EcoRI");
  assert.equal(state.rightEnzyme, "BamHI");
  assert.deepEqual(state.selectedEnzymes, ["EcoRI", "BamHI"]);
  assert.equal(state.vectorFragmentIndex, "2");
  assert.equal(state.insertFragmentIndex, "1");
});

test("enzyme map clicks update selection state without side effects", () => {
  const firstClick = selection.applyEnzymeClick(selection.createState({}), "EcoRI");
  assert.equal(firstClick.state.leftEnzyme, "EcoRI");
  assert.equal(firstClick.state.rightEnzyme, "EcoRI");

  const secondClick = selection.applyEnzymeClick(firstClick.state, "BamHI");
  assert.equal(secondClick.state.leftEnzyme, "EcoRI");
  assert.equal(secondClick.state.rightEnzyme, "BamHI");
  assert.match(secondClick.message, /locally/i);

  const clearRight = selection.applyEnzymeClick(
    selection.createState({ leftEnzyme: "EcoRI", rightEnzyme: "EcoRI" }),
    "EcoRI"
  );
  assert.equal(clearRight.state.leftEnzyme, "EcoRI");
  assert.equal(clearRight.state.rightEnzyme, "");
});

test("fragment toggles accept Django field names", () => {
  const initial = selection.createState({ vectorFragmentIndex: "1" });
  const cleared = selection.toggleFragment(initial, "vector_fragment_index", "1");
  assert.equal(cleared.vectorFragmentIndex, "");

  const selected = selection.toggleFragment(cleared, "insert_fragment_index", "3");
  assert.equal(selected.insertFragmentIndex, "3");
});

test("authoritative selection key includes selected enzyme changes", () => {
  const base = selection.createState({
    vectorAsset: "sequence_file:1:vec1",
    insertAsset: "sequence_file:2:ins1",
    assemblyStrategy: "restriction_ligation",
    leftEnzyme: "EcoRI",
    rightEnzyme: "BamHI",
  });
  const overlayOnly = selection.setSelectedEnzymes(base, ["EcoRI", "BamHI", "BbsI"]);

  assert.notEqual(selection.selectionKey(base), selection.selectionKey(overlayOnly));
  assert.notEqual(
    selection.selectionKey(base),
    selection.selectionKey(selection.setField(base, "right_enzyme", "BbsI"))
  );
});

test("selected enzyme changes derive backend enzyme pair and clear stale fragments", () => {
  const state = selection.setSelectedEnzymes(
    selection.createState({
      leftEnzyme: "EcoRI",
      rightEnzyme: "EcoRI",
      vectorFragmentIndex: "2",
      insertFragmentIndex: "1",
    }),
    ["BsaI"]
  );

  assert.equal(state.leftEnzyme, "BsaI");
  assert.equal(state.rightEnzyme, "BsaI");
  assert.deepEqual(state.selectedEnzymes, ["BsaI"]);
  assert.equal(state.vectorFragmentIndex, "");
  assert.equal(state.insertFragmentIndex, "");
});

test("preview interaction script does not programmatically submit on selection changes", () => {
  const source = fs.readFileSync(
    path.join(__dirname, "cloning_assembly_preview.js"),
    "utf8"
  );

  assert.doesNotMatch(source, /requestSubmit/);
  assert.doesNotMatch(source, /HTMLFormElement\.prototype\.submit/);
  assert.doesNotMatch(source, /submitPreviewWhenReady/);
});

test("preview interaction script paginates the enzyme selection table", () => {
  const source = fs.readFileSync(
    path.join(__dirname, "cloning_assembly_preview.js"),
    "utf8"
  );

  assert.match(source, /ENZYME_TABLE_PAGE_SIZE/);
  assert.match(source, /data-cloning-enzyme-prev-page/);
  assert.match(source, /data-cloning-enzyme-next-page/);
  assert.match(source, /data-cloning-enzyme-page-info/);
});
