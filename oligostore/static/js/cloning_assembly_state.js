(function (root, factory) {
  if (typeof module === "object" && module.exports) {
    module.exports = factory();
  } else {
    root.RestrictionLigationSelection = factory();
  }
})(typeof globalThis !== "undefined" ? globalThis : this, function () {
  function normalizeText(value) {
    return String(value || "").trim();
  }

  function normalizeList(values) {
    const output = [];
    (Array.isArray(values) ? values : []).forEach((value) => {
      const normalized = normalizeText(value);
      if (normalized && !output.includes(normalized)) {
        output.push(normalized);
      }
    });
    return output;
  }

  function normalizeFieldName(fieldName) {
    const normalized = normalizeText(fieldName);
    if (normalized === "vector_fragment_index") {
      return "vectorFragmentIndex";
    }
    if (normalized === "insert_fragment_index") {
      return "insertFragmentIndex";
    }
    if (normalized === "left_enzyme") {
      return "leftEnzyme";
    }
    if (normalized === "right_enzyme") {
      return "rightEnzyme";
    }
    if (normalized === "is_circular") {
      return "isCircular";
    }
    return normalized;
  }

  function createState(initial) {
    const source = initial || {};
    const leftEnzyme = normalizeText(source.leftEnzyme || source.left_enzyme);
    const rightEnzyme = normalizeText(source.rightEnzyme || source.right_enzyme);
    const selectedEnzymes = normalizeList(source.selectedEnzymes || source.selected_enzymes);
    return {
      vectorAsset: normalizeText(source.vectorAsset || source.vector_asset),
      insertAsset: normalizeText(source.insertAsset || source.insert_asset),
      assemblyStrategy: normalizeText(source.assemblyStrategy || source.assembly_strategy),
      isCircular: normalizeText(source.isCircular || source.is_circular),
      leftEnzyme,
      rightEnzyme,
      selectedEnzymes: selectedEnzymes.length ? selectedEnzymes : normalizeList([leftEnzyme, rightEnzyme]),
      vectorFragmentIndex: normalizeText(source.vectorFragmentIndex || source.vector_fragment_index),
      insertFragmentIndex: normalizeText(source.insertFragmentIndex || source.insert_fragment_index),
    };
  }

  function setField(state, fieldName, value) {
    const normalizedField = normalizeFieldName(fieldName);
    return createState({
      ...(state || {}),
      [normalizedField]: value,
    });
  }

  function setSelectedEnzymes(state, values) {
    const selectedEnzymes = normalizeList(values);
    return createState({
      ...(state || {}),
      selectedEnzymes,
      leftEnzyme: selectedEnzymes[0] || "",
      rightEnzyme: selectedEnzymes[1] || selectedEnzymes[0] || "",
      vectorFragmentIndex: "",
      insertFragmentIndex: "",
    });
  }

  function applyEnzymeClick(state, enzymeName) {
    const enzyme = normalizeText(enzymeName);
    const current = createState(state);
    if (!enzyme) {
      return {
        state: current,
        message: "That restriction site is not available for the selected vector.",
        tone: "text-error",
      };
    }

    let left = current.leftEnzyme;
    let right = current.rightEnzyme;
    let message = "Selection preview updated locally. Click Preview Assembly to run authoritative validation.";

    if (left === enzyme && right === enzyme) {
      right = "";
      message = "Right enzyme cleared locally. Choose another site or preview the current pair.";
    } else if (!left) {
      left = enzyme;
      if (!right) {
        right = enzyme;
      }
    } else if (left === enzyme && right !== enzyme) {
      right = enzyme;
    } else if (right === enzyme && left !== enzyme) {
      left = enzyme;
    } else if (left === right) {
      right = enzyme;
    } else {
      right = enzyme;
    }

    return {
      state: createState({
        ...current,
        leftEnzyme: left,
        rightEnzyme: right,
      }),
      message,
      tone: "text-warning",
    };
  }

  function toggleFragment(state, fieldName, fragmentIndex) {
    const normalizedField = normalizeFieldName(fieldName);
    const normalizedIndex = normalizeText(fragmentIndex);
    const current = createState(state);
    if (!normalizedField || !Object.prototype.hasOwnProperty.call(current, normalizedField)) {
      return current;
    }
    return setField(
      current,
      normalizedField,
      current[normalizedField] === normalizedIndex ? "" : normalizedIndex
    );
  }

  function selectionKey(state) {
    const current = createState(state);
    return [
      current.vectorAsset,
      current.insertAsset,
      current.assemblyStrategy,
      current.isCircular,
      current.leftEnzyme,
      current.rightEnzyme,
      current.selectedEnzymes.join(","),
      current.vectorFragmentIndex,
      current.insertFragmentIndex,
    ].join("|");
  }

  return {
    applyEnzymeClick,
    createState,
    selectionKey,
    setField,
    setSelectedEnzymes,
    toggleFragment,
  };
});
