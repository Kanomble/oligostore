(function () {
  const fragmentAnchors = new Map();
  const selectionStore = window.RestrictionLigationSelection;
  let selectionState = null;
  let authoritativeSelectionKey = "";
  let pendingRenderFrame = null;
  let enzymeTablePage = 1;
  const ENZYME_TABLE_PAGE_SIZE = 12;
  const SVG_NS = "http://www.w3.org/2000/svg";

  function readMapPayload() {
    const script = document.getElementById("cloning-map-data");
    if (!script) {
      return null;
    }
    try {
      return JSON.parse(script.textContent || "{}");
    } catch (_error) {
      return null;
    }
  }

  const mapPayload = readMapPayload();

  function getCloningPreviewForm() {
    return document.getElementById("cloning-assembly-preview-form");
  }

  function getCloningSaveForm() {
    return document.querySelector("[data-cloning-save-form]");
  }

  function setCloningVisualFeedback(message, tone) {
    const feedback = document.querySelector("[data-cloning-visual-feedback]");
    if (!feedback) {
      return;
    }
    feedback.textContent = message || "";
    feedback.classList.remove("text-error", "text-success", "text-warning");
    if (tone) {
      feedback.classList.add(tone);
    }
  }

  function setFormLoading(form, message) {
    if (!form) {
      return;
    }
    form.classList.add("is-loading");
    form.setAttribute("aria-busy", "true");
    form.querySelectorAll('button[type="submit"]').forEach((button) => {
      button.classList.add("loading");
      button.disabled = true;
    });
    if (message) {
      setCloningVisualFeedback(message, "text-warning");
    }
  }

  function selectedValues(select) {
    if (!select || !select.options) {
      return [];
    }
    return Array.from(select.options)
      .filter((option) => option.selected && option.value)
      .map((option) => option.value);
  }

  function readSelectionState(form) {
    if (!selectionStore || !form) {
      return null;
    }
    return selectionStore.createState({
      vectorAsset: form.elements.vector_asset ? form.elements.vector_asset.value : "",
      insertAsset: form.elements.insert_asset ? form.elements.insert_asset.value : "",
      assemblyStrategy: form.elements.assembly_strategy ? form.elements.assembly_strategy.value : "",
      isCircular: form.elements.is_circular ? form.elements.is_circular.value : "",
      leftEnzyme: form.elements.left_enzyme ? form.elements.left_enzyme.value : "",
      rightEnzyme: form.elements.right_enzyme ? form.elements.right_enzyme.value : "",
      selectedEnzymes: selectedValues(form.elements.selected_enzymes),
      vectorFragmentIndex: form.elements.vector_fragment_index ? form.elements.vector_fragment_index.value : "",
      insertFragmentIndex: form.elements.insert_fragment_index ? form.elements.insert_fragment_index.value : "",
    });
  }

  function setMaybeField(form, name, value) {
    const field = form && form.elements[name];
    if (field) {
      field.value = value || "";
    }
  }

  function writeSelectionStateToForm(form, state) {
    if (!form || !state) {
      return;
    }
    setMaybeField(form, "left_enzyme", state.leftEnzyme);
    setMaybeField(form, "right_enzyme", state.rightEnzyme);
    setMaybeField(form, "is_circular", state.isCircular);
    setMaybeField(form, "vector_fragment_index", state.vectorFragmentIndex);
    setMaybeField(form, "insert_fragment_index", state.insertFragmentIndex);
    const selected = form.elements.selected_enzymes;
    if (selected && selected.options) {
      const selectedSet = new Set(state.selectedEnzymes);
      Array.from(selected.options).forEach((option) => {
        option.selected = selectedSet.has(option.value);
      });
    }
  }

  function writeSelectionStateToSaveForm(state) {
    const saveForm = getCloningSaveForm();
    if (!saveForm || !state) {
      return;
    }
    setMaybeField(saveForm, "left_enzyme", state.leftEnzyme);
    setMaybeField(saveForm, "right_enzyme", state.rightEnzyme);
    setMaybeField(saveForm, "is_circular", state.isCircular);
    setMaybeField(saveForm, "vector_fragment_index", state.vectorFragmentIndex);
    setMaybeField(saveForm, "insert_fragment_index", state.insertFragmentIndex);
  }

  function refreshVisualEnzymeSelection(form) {
    if (!form) {
      return;
    }
    const selectedSet = new Set(selectedValues(form.elements.selected_enzymes));
    document.querySelectorAll("[data-cloning-enzyme-picker]").forEach((site) => {
      const enzymeName = site.dataset.enzymeName || "";
      site.classList.toggle("is-selected", Boolean(enzymeName && selectedSet.has(enzymeName)));
    });
  }

  function updateStrategySummary(state) {
    const target = document.querySelector("[data-cloning-selected-fragment-summary]");
    if (!target) {
      return;
    }
    const vector = state && state.vectorFragmentIndex ? `Vector fragment ${state.vectorFragmentIndex}` : "No vector fragment";
    const insert = state && state.insertFragmentIndex ? `Insert fragment ${state.insertFragmentIndex}` : "no insert fragment";
    target.textContent = `${vector}; ${insert} selected.`;
  }

  function fragmentFieldForRole(role) {
    if (role === "vector") {
      return "vectorFragmentIndex";
    }
    if (role === "insert") {
      return "insertFragmentIndex";
    }
    return "";
  }

  function refreshVisualFragmentSelection(state) {
    if (!state) {
      return;
    }
    document.querySelectorAll("[data-cloning-region-picker], [data-cloning-digest-fragment]").forEach((fragment) => {
      const stateField = fragmentFieldForRole(fragment.dataset.mapRole || "");
      const selectedValue = stateField ? state[stateField] : "";
      const isSelected = Boolean(selectedValue && selectedValue === (fragment.dataset.fragmentIndex || ""));
      fragment.classList.toggle("is-selected", isSelected);
      fragment.setAttribute("data-region-selected", isSelected ? "1" : "0");
    });
  }

  function updateSaveAvailability(state) {
    const saveForm = getCloningSaveForm();
    if (!saveForm || !selectionStore || !state) {
      return;
    }
    const saveButton = saveForm.querySelector('button[type="submit"]');
    const isStale = selectionStore.selectionKey(state) !== authoritativeSelectionKey;
    if (saveButton) {
      saveButton.disabled = isStale;
      saveButton.classList.toggle("btn-disabled", isStale);
      saveButton.title = isStale ? "Run Preview Assembly again before saving this changed selection." : "";
    }
    saveForm.dataset.previewStale = isStale ? "1" : "0";
    document.querySelectorAll("[data-cloning-authoritative-preview]").forEach((preview) => {
      preview.hidden = isStale;
    });
  }

  function renderLocalSelection(form, message, tone) {
    if (!form || !selectionState) {
      return;
    }
    writeSelectionStateToForm(form, selectionState);
    writeSelectionStateToSaveForm(selectionState);
    refreshVisualEnzymeSelection(form);
    refreshVisualFragmentSelection(selectionState);
    updateStrategySummary(selectionState);
    renderDynamicMaps(form, selectionState);
    refreshUnifiedEnzymeTable(form);
    updateSaveAvailability(selectionState);
    if (message) {
      setCloningVisualFeedback(message, tone || "text-warning");
    }
  }

  function scheduleLocalSelectionRender(form, message, tone) {
    if (pendingRenderFrame !== null && typeof window.cancelAnimationFrame === "function") {
      window.cancelAnimationFrame(pendingRenderFrame);
    }
    if (typeof window.requestAnimationFrame !== "function") {
      renderLocalSelection(form, message, tone);
      return;
    }
    pendingRenderFrame = window.requestAnimationFrame(() => {
      pendingRenderFrame = null;
      renderLocalSelection(form, message, tone);
    });
  }

  function applySelectionUpdate(form, nextState, message, tone) {
    if (!selectionStore || !form || !nextState) {
      return;
    }
    selectionState = selectionStore.createState(nextState);
    scheduleLocalSelectionRender(form, message, tone);
  }

  function clearFragmentAnchors(role) {
    const roles = role ? [role] : Array.from(fragmentAnchors.keys());
    roles.forEach((mapRole) => {
      const previous = fragmentAnchors.get(mapRole);
      if (previous && previous.element) {
        previous.element.classList.remove("is-fragment-anchor");
      }
      fragmentAnchors.delete(mapRole);
    });
  }

  function siteData(trigger) {
    const map = trigger.closest("[data-cloning-sequence-map]");
    const position = Number.parseInt(trigger.dataset.sitePosition || "", 10);
    return {
      element: trigger,
      role: trigger.dataset.mapRole || (map ? map.dataset.mapRole : ""),
      mapShape: map ? map.dataset.mapShape || "linear" : "linear",
      enzymeName: trigger.dataset.enzymeName || "",
      siteId: trigger.dataset.siteId || "",
      position: Number.isFinite(position) ? position : null,
    };
  }

  function fragmentData(trigger) {
    const start = Number.parseInt(trigger.dataset.fragmentStart || "", 10);
    const end = Number.parseInt(trigger.dataset.fragmentEnd || "", 10);
    return {
      element: trigger,
      role: trigger.dataset.mapRole || "",
      fieldName: trigger.dataset.fragmentField || "",
      index: trigger.dataset.fragmentIndex || "",
      start: Number.isFinite(start) ? start : null,
      end: Number.isFinite(end) ? end : null,
      wrapsOrigin: trigger.dataset.fragmentWrapsOrigin === "1",
    };
  }

  function enzymePayloadByName() {
    const byName = new Map();
    ((mapPayload && mapPayload.enzymes) || []).forEach((enzyme) => {
      if (enzyme && enzyme.name) {
        byName.set(enzyme.name, enzyme);
      }
    });
    return byName;
  }

  function selectedEnzymePayloads(state) {
    const byName = enzymePayloadByName();
    return (state && state.selectedEnzymes ? state.selectedEnzymes : [])
      .map((name) => byName.get(name))
      .filter(Boolean);
  }

  function cutPositionsForRole(state, role) {
    const fieldName = role === "vector" ? "vectorCutPositions" : "insertCutPositions";
    const cuts = [];
    selectedEnzymePayloads(state).forEach((enzyme) => {
      (enzyme[fieldName] || []).forEach((position) => {
        const numericPosition = Number.parseInt(position, 10);
        if (Number.isFinite(numericPosition)) {
          cuts.push({
            enzymeName: enzyme.name,
            siteSequence: enzyme.siteSequence || "",
            position: numericPosition,
          });
        }
      });
    });
    return cuts.sort((a, b) => a.position - b.position || a.enzymeName.localeCompare(b.enzymeName));
  }

  function digestFragmentsForRole(state, role) {
    const sequenceMeta = mapPayload ? mapPayload[role] : null;
    const sequenceLength = sequenceMeta ? Number.parseInt(sequenceMeta.sequenceLength, 10) : 0;
    const cuts = cutPositionsForRole(state, role);
    const uniqueCuts = Array.from(new Set(cuts.map((cut) => cut.position))).sort((a, b) => a - b);
    if (!sequenceLength || uniqueCuts.length === 0) {
      return [];
    }
    if (sequenceMeta.mapShape === "circular") {
      return uniqueCuts.map((start, index) => {
        const end = uniqueCuts[(index + 1) % uniqueCuts.length];
        return {
          index: index + 1,
          start,
          end,
          wrapsOrigin: end <= start,
          length: end <= start ? sequenceLength - start + end : end - start,
        };
      });
    }
    const boundaries = [0].concat(uniqueCuts, [sequenceLength]);
    const fragments = [];
    for (let index = 0; index < boundaries.length - 1; index += 1) {
      const start = boundaries[index];
      const end = boundaries[index + 1];
      if (end > start) {
        fragments.push({
          index: fragments.length + 1,
          start,
          end,
          wrapsOrigin: false,
          length: end - start,
        });
      }
    }
    return fragments;
  }

  function roleFieldName(role) {
    return role === "vector" ? "vector_fragment_index" : "insert_fragment_index";
  }

  function createSvgElement(tagName, attributes) {
    const element = document.createElementNS(SVG_NS, tagName);
    Object.entries(attributes || {}).forEach(([name, value]) => {
      element.setAttribute(name, String(value));
    });
    return element;
  }

  function polarPoint(position, sequenceLength, radius) {
    const angle = sequenceLength > 0 ? ((position / sequenceLength) * Math.PI * 2) - (Math.PI / 2) : -Math.PI / 2;
    return {
      x: (50 + radius * Math.cos(angle)).toFixed(3),
      y: (50 + radius * Math.sin(angle)).toFixed(3),
      anchor: Math.cos(angle) > 0.25 ? "start" : Math.cos(angle) < -0.25 ? "end" : "middle",
    };
  }

  function fragmentTitle(fragment) {
    if (fragment.wrapsOrigin) {
      return `Fragment ${fragment.index}: ${fragment.start + 1}-${fragment.sequenceLength || ""} + 1-${fragment.end}, ${fragment.length} bp`;
    }
    return `Fragment ${fragment.index}: ${fragment.start + 1}-${fragment.end}, ${fragment.length} bp`;
  }

  function applyFragmentDataset(element, role, fragment) {
    element.dataset.cloningRegionPicker = "";
    element.dataset.mapRole = role;
    element.dataset.fragmentField = roleFieldName(role);
    element.dataset.fragmentIndex = String(fragment.index);
    element.dataset.fragmentStart = String(fragment.start);
    element.dataset.fragmentEnd = String(fragment.end);
    element.dataset.fragmentWrapsOrigin = fragment.wrapsOrigin ? "1" : "0";
    element.dataset.regionSelected = "0";
  }

  function renderCircularRegion(svg, role, fragment, sequenceLength) {
    const segments = fragment.wrapsOrigin ? [[fragment.start, sequenceLength], [0, fragment.end]] : [[fragment.start, fragment.end]];
    segments.forEach(([start, end]) => {
      if (end <= start) {
        return;
      }
      const segmentLength = Math.max(1, ((end - start) / sequenceLength) * 100);
      const region = createSvgElement("circle", {
        class: "assembly-plasmid-region is-selectable",
        cx: 50,
        cy: 50,
        r: 35,
        pathLength: 100,
        "stroke-dasharray": `${segmentLength.toFixed(3)} 100`,
        "stroke-dashoffset": `${(-(start / sequenceLength) * 100).toFixed(3)}`,
        transform: "rotate(-90 50 50)",
        tabindex: 0,
        role: "button",
      });
      fragment.sequenceLength = sequenceLength;
      applyFragmentDataset(region, role, fragment);
      const title = createSvgElement("title");
      title.textContent = fragmentTitle(fragment);
      region.appendChild(title);
      svg.appendChild(region);
    });
  }

  function renderCircularMap(panel, role, state) {
    const meta = mapPayload[role];
    const sequenceLength = Number.parseInt(meta.sequenceLength, 10);
    const body = panel.querySelector("[data-cloning-map-body]");
    if (!body) {
      return;
    }
    body.replaceChildren();
    const wrap = document.createElement("div");
    wrap.className = "assembly-plasmid-wrap";
    const svg = createSvgElement("svg", {
      class: "assembly-plasmid-map",
      viewBox: "-8 -8 116 116",
      role: "img",
      "aria-label": `Circular ${role} map`,
    });
    svg.appendChild(createSvgElement("circle", { class: "assembly-plasmid-backbone", cx: 50, cy: 50, r: 35 }));
    digestFragmentsForRole(state, role).forEach((fragment) => renderCircularRegion(svg, role, fragment, sequenceLength));
    cutPositionsForRole(state, role).forEach((cut) => {
      const markerPoint = polarPoint(cut.position, sequenceLength, 38);
      const labelPoint = polarPoint(cut.position, sequenceLength, 47);
      const marker = createSvgElement("g", {
        class: "assembly-site-marker is-selected",
        tabindex: 0,
        role: "button",
        "data-cloning-enzyme-picker": "",
        "data-site-id": `${role}-${cut.enzymeName.toLowerCase()}-${cut.position}`,
        "data-map-role": role,
        "data-enzyme-name": cut.enzymeName,
        "data-enzyme-title": `${cut.enzymeName} at ${role} base ${cut.position + 1}`,
        "data-site-position": cut.position,
        "data-sequence-length": sequenceLength,
      });
      const title = createSvgElement("title");
      title.textContent = `${cut.enzymeName} at ${role} base ${cut.position + 1}`;
      marker.appendChild(title);
      marker.appendChild(createSvgElement("line", { x1: 50, y1: 50, x2: markerPoint.x, y2: markerPoint.y }));
      marker.appendChild(createSvgElement("circle", { cx: markerPoint.x, cy: markerPoint.y, r: 2.7 }));
      const text = createSvgElement("text", { x: labelPoint.x, y: labelPoint.y, "text-anchor": labelPoint.anchor });
      text.textContent = cut.enzymeName;
      marker.appendChild(text);
      svg.appendChild(marker);
    });
    wrap.appendChild(svg);
    body.appendChild(wrap);
    if (!state.selectedEnzymes.length) {
      appendMapStatus(body, "Select enzymes to show digest fragments.");
    } else if (!cutPositionsForRole(state, role).length) {
      appendMapStatus(body, "Selected enzymes produce no cut sites on this sequence.");
    }
  }

  function appendMapStatus(body, message) {
    const status = document.createElement("div");
    status.className = "assembly-map-status";
    status.textContent = message;
    body.appendChild(status);
  }

  function renderLinearMap(panel, role, state) {
    const meta = mapPayload[role];
    const sequenceLength = Number.parseInt(meta.sequenceLength, 10);
    const body = panel.querySelector("[data-cloning-map-body]");
    if (!body) {
      return;
    }
    body.replaceChildren();
    const map = document.createElement("div");
    map.className = "assembly-linear-map";
    map.setAttribute("aria-label", `Linear ${role} map`);
    const track = document.createElement("div");
    track.className = "assembly-linear-track";
    digestFragmentsForRole(state, role).forEach((fragment) => {
      const region = document.createElement("button");
      region.type = "button";
      region.className = "assembly-linear-region is-selectable";
      region.style.left = `${((fragment.start / sequenceLength) * 100).toFixed(3)}%`;
      region.style.width = `${Math.max(1.25, (fragment.length / sequenceLength) * 100).toFixed(3)}%`;
      region.title = fragmentTitle(fragment);
      applyFragmentDataset(region, role, fragment);
      track.appendChild(region);
    });
    cutPositionsForRole(state, role).forEach((cut) => {
      const site = document.createElement("button");
      site.type = "button";
      site.className = "assembly-linear-site is-selected";
      site.style.left = `${((cut.position / sequenceLength) * 100).toFixed(3)}%`;
      site.title = `${cut.enzymeName} at ${role} base ${cut.position + 1}`;
      site.dataset.cloningEnzymePicker = "";
      site.dataset.siteId = `${role}-${cut.enzymeName.toLowerCase()}-${cut.position}`;
      site.dataset.mapRole = role;
      site.dataset.enzymeName = cut.enzymeName;
      site.dataset.enzymeTitle = site.title;
      site.dataset.sitePosition = String(cut.position);
      site.dataset.sequenceLength = String(sequenceLength);
      const label = document.createElement("span");
      label.className = "assembly-linear-site-label";
      label.textContent = cut.enzymeName;
      site.appendChild(label);
      track.appendChild(site);
    });
    map.appendChild(track);
    body.appendChild(map);
    if (!state.selectedEnzymes.length) {
      appendMapStatus(body, "Select enzymes to show digest fragments.");
    } else if (!cutPositionsForRole(state, role).length) {
      appendMapStatus(body, "Selected enzymes produce no cut sites on this sequence.");
    }
  }

  function renderDigestFragments(role, state) {
    const container = document.querySelector(`[data-cloning-digest-fragments][data-map-role="${role}"]`);
    if (!container || !mapPayload || !mapPayload[role]) {
      return;
    }
    const fragments = digestFragmentsForRole(state, role);
    container.replaceChildren();
    if (!state.selectedEnzymes.length) {
      const empty = document.createElement("div");
      empty.className = "text-sm opacity-70";
      empty.textContent = "Select enzymes to preview digest fragments.";
      container.appendChild(empty);
      return;
    }
    if (!fragments.length) {
      const empty = document.createElement("div");
      empty.className = "text-sm opacity-70";
      empty.textContent = "Selected enzymes produce no cut sites on this sequence.";
      container.appendChild(empty);
      return;
    }
    const wrap = document.createElement("div");
    wrap.className = "overflow-x-auto";
    const table = document.createElement("table");
    table.className = "assembly-digest-table";
    table.innerHTML = "<thead><tr><th>Fragment</th><th>Coordinates</th><th>Length</th><th>Boundaries</th></tr></thead>";
    const tbody = document.createElement("tbody");
    fragments.forEach((fragment) => {
      const row = document.createElement("tr");
      row.dataset.cloningDigestFragment = "";
      row.dataset.mapRole = role;
      row.dataset.fragmentField = roleFieldName(role);
      row.dataset.fragmentIndex = String(fragment.index);
      row.dataset.fragmentStart = String(fragment.start);
      row.dataset.fragmentEnd = String(fragment.end);
      row.dataset.fragmentWrapsOrigin = fragment.wrapsOrigin ? "1" : "0";
      const coordinates = fragment.wrapsOrigin
        ? `${fragment.start + 1}-${mapPayload[role].sequenceLength} + 1-${fragment.end}`
        : `${fragment.start + 1}-${fragment.end}`;
      row.innerHTML = `<td>${fragment.index}</td><td>${coordinates}</td><td>${fragment.length} bp</td><td>selected cuts</td>`;
      tbody.appendChild(row);
    });
    table.appendChild(tbody);
    wrap.appendChild(table);
    container.appendChild(wrap);
  }

  function renderDynamicMaps(form, state) {
    if (!mapPayload || !state) {
      return;
    }
    ["vector", "insert"].forEach((role) => {
      const panel = document.querySelector(`[data-cloning-sequence-map][data-map-role="${role}"]`);
      if (!panel || !mapPayload[role]) {
        return;
      }
      if (mapPayload[role].mapShape === "circular") {
        renderCircularMap(panel, role, state);
      } else {
        renderLinearMap(panel, role, state);
      }
      renderDigestFragments(role, state);
    });
    refreshVisualFragmentSelection(state);
    refreshVisualEnzymeSelection(form);
  }

  function refreshUnifiedEnzymeTable(form) {
    const selectedSet = new Set(selectedValues(form && form.elements.selected_enzymes));
    const byName = enzymePayloadByName();
    const rows = Array.from(document.querySelectorAll("[data-cloning-enzyme-row]"));
    const searchInput = document.querySelector("[data-cloning-enzyme-search]");
    const query = String(searchInput ? searchInput.value || "" : "").trim().toLowerCase();
    const matchingRows = rows.filter((row) => {
      const name = String(row.dataset.enzymeName || "").toLowerCase();
      return !query || name.includes(query);
    });
    const totalPages = Math.max(1, Math.ceil(matchingRows.length / ENZYME_TABLE_PAGE_SIZE));
    enzymeTablePage = Math.min(Math.max(1, enzymeTablePage), totalPages);
    const pageStart = (enzymeTablePage - 1) * ENZYME_TABLE_PAGE_SIZE;
    const visibleRows = new Set(matchingRows.slice(pageStart, pageStart + ENZYME_TABLE_PAGE_SIZE));

    rows.forEach((row) => {
      const name = row.dataset.enzymeName || "";
      const checkbox = row.querySelector("[data-cloning-enzyme-toggle]");
      const enzyme = byName.get(name);
      row.classList.toggle("is-selected", selectedSet.has(name));
      row.hidden = !visibleRows.has(row);
      if (checkbox) {
        checkbox.checked = selectedSet.has(name);
      }
      const vectorTarget = row.querySelector("[data-cloning-enzyme-vector-count]");
      const insertTarget = row.querySelector("[data-cloning-enzyme-insert-count]");
      if (vectorTarget) {
        vectorTarget.textContent = enzyme ? String((enzyme.vectorCutPositions || []).length) : "0";
      }
      if (insertTarget) {
        insertTarget.textContent = enzyme ? String((enzyme.insertCutPositions || []).length) : "0";
      }
    });
    const summary = document.querySelector("[data-cloning-enzyme-selection-summary]");
    if (summary) {
      const selected = Array.from(selectedSet);
      const visibleEnd = Math.min(pageStart + ENZYME_TABLE_PAGE_SIZE, matchingRows.length);
      const rangeLabel = matchingRows.length
        ? `Showing ${pageStart + 1}-${visibleEnd} of ${matchingRows.length} enzymes.`
        : "No matching enzymes.";
      summary.textContent = selected.length
        ? `${rangeLabel} Selected: ${selected.join(", ")}`
        : `${rangeLabel} No enzymes selected.`;
    }
    const pageInfo = document.querySelector("[data-cloning-enzyme-page-info]");
    if (pageInfo) {
      pageInfo.textContent = `Page ${enzymeTablePage} of ${totalPages}`;
    }
    const prevButton = document.querySelector("[data-cloning-enzyme-prev-page]");
    if (prevButton) {
      prevButton.disabled = enzymeTablePage <= 1;
      prevButton.classList.toggle("btn-disabled", enzymeTablePage <= 1);
    }
    const nextButton = document.querySelector("[data-cloning-enzyme-next-page]");
    if (nextButton) {
      nextButton.disabled = enzymeTablePage >= totalPages;
      nextButton.classList.toggle("btn-disabled", enzymeTablePage >= totalPages);
    }
  }

  function setSelectedEnzymeValues(form, values) {
    const select = form && form.elements.selected_enzymes;
    if (!select || !select.options) {
      return;
    }
    const selectedSet = new Set(values);
    Array.from(select.options).forEach((option) => {
      option.selected = selectedSet.has(option.value);
    });
    select.dispatchEvent(new Event("change", { bubbles: true }));
  }

  function fragmentMatchesPair(fragment, firstSite, secondSite) {
    if (fragment.start === null || fragment.end === null || firstSite.position === null || secondSite.position === null) {
      return false;
    }
    if (firstSite.mapShape === "circular") {
      const forwardMatch = (
        fragment.start === firstSite.position &&
        fragment.end === secondSite.position &&
        fragment.wrapsOrigin === (secondSite.position <= firstSite.position)
      );
      const reverseMatch = (
        fragment.start === secondSite.position &&
        fragment.end === firstSite.position &&
        fragment.wrapsOrigin === (firstSite.position <= secondSite.position)
      );
      return forwardMatch || reverseMatch;
    }
    return (
      fragment.start === Math.min(firstSite.position, secondSite.position) &&
      fragment.end === Math.max(firstSite.position, secondSite.position)
    );
  }

  function matchingFragmentForPair(form, firstSite, secondSite) {
    const regionCandidates = Array.from(
      document.querySelectorAll("[data-cloning-region-picker]")
    )
      .filter((candidate) => candidate.dataset.mapRole === firstSite.role)
      .map(fragmentData);
    const regionMatch = regionCandidates.find((fragment) => fragmentMatchesPair(fragment, firstSite, secondSite));
    if (regionMatch) {
      return regionMatch;
    }

    return Array.from(
      document.querySelectorAll("[data-cloning-digest-fragment]")
    )
      .filter((candidate) => candidate.dataset.mapRole === firstSite.role)
      .map(fragmentData)
      .find((fragment) => fragmentMatchesPair(fragment, firstSite, secondSite)) || null;
  }

  function handleShiftClickSite(trigger, form) {
    const current = siteData(trigger);
    if (!current.role || !current.enzymeName || current.position === null) {
      setCloningVisualFeedback("That restriction site is missing fragment-selection metadata.", "text-error");
      return;
    }

    const previous = fragmentAnchors.get(current.role);
    if (!previous || previous.siteId === current.siteId) {
      clearFragmentAnchors(current.role);
      fragmentAnchors.set(current.role, current);
      trigger.classList.add("is-fragment-anchor");
      setCloningVisualFeedback(`Fragment anchor set on ${current.role} ${current.enzymeName}. Shift-click the second site on the same map.`, "text-warning");
      return;
    }

    const fragment = matchingFragmentForPair(form, previous, current);
    if (!fragment || !fragment.fieldName || !fragment.index) {
      clearFragmentAnchors(current.role);
      setCloningVisualFeedback("No selectable digest fragment matches those two sites. Use adjacent same-enzyme cut sites or select the fragment from the map.", "text-error");
      return;
    }

    const field = form.elements[fragment.fieldName];
    if (!field) {
      clearFragmentAnchors(current.role);
      setCloningVisualFeedback("That fragment cannot be submitted with the current enzyme choices.", "text-error");
      return;
    }

    if (selectionStore) {
      let nextState = selectionState || readSelectionState(form);
      nextState = selectionStore.setField(
        nextState,
        fragment.fieldName === "vector_fragment_index" ? "vectorFragmentIndex" : "insertFragmentIndex",
        fragment.index
      );
      applySelectionUpdate(
        form,
        nextState,
        `${current.role} fragment ${fragment.index} selected locally. Click Preview Assembly to validate.`,
        "text-warning"
      );
    } else {
      field.value = fragment.index;
      refreshVisualEnzymeSelection(form);
      refreshVisualFragmentSelection(readSelectionState(form));
    }
    clearFragmentAnchors(current.role);
  }

  document.addEventListener("click", function (event) {
    const regionTrigger = event.target.closest("[data-cloning-region-picker], [data-cloning-fragment-picker]");
    if (regionTrigger) {
      const form = getCloningPreviewForm();
      if (!form) {
        return;
      }
      const field = form.elements[regionTrigger.dataset.fragmentField];
      if (!field) {
        return;
      }
      const nextValue = regionTrigger.dataset.fragmentIndex || "";
      if (selectionStore) {
        applySelectionUpdate(
          form,
          selectionStore.toggleFragment(selectionState || readSelectionState(form), regionTrigger.dataset.fragmentField, nextValue),
          "Fragment selection updated locally. Click Preview Assembly to validate the ligation.",
          "text-warning"
        );
      } else {
        field.value = field.value === nextValue ? "" : nextValue;
        refreshVisualFragmentSelection(readSelectionState(form));
      }
      return;
    }

    const enzymeTrigger = event.target.closest("[data-cloning-enzyme-picker]");
    if (!enzymeTrigger) {
      return;
    }

    const form = getCloningPreviewForm();
    if (!form) {
      return;
    }

    if (event.shiftKey) {
      handleShiftClickSite(enzymeTrigger, form);
      return;
    }

    const current = siteData(enzymeTrigger);
    if (!current.role || !current.enzymeName || current.position === null) {
      setCloningVisualFeedback("That restriction site is missing fragment-selection metadata.", "text-error");
      return;
    }
    const previous = fragmentAnchors.get(current.role);
    if (previous && previous.siteId === current.siteId) {
      clearFragmentAnchors(current.role);
      setCloningVisualFeedback("Cut-site anchor cleared.", "text-warning");
      return;
    }
    clearFragmentAnchors();
    clearFragmentAnchors(current.role);
    fragmentAnchors.set(current.role, current);
    enzymeTrigger.classList.add("is-fragment-anchor");
    setCloningVisualFeedback(`Cut-site anchor set on ${current.role}. Shift-click another cut site or click a fragment.`, "text-warning");
  });

  document.addEventListener("change", function (event) {
    const form = getCloningPreviewForm();
    if (!form || !selectionStore || !event.target) {
      return;
    }
    const target = event.target;
    if (target.name === "left_enzyme") {
      applySelectionUpdate(
        form,
        selectionStore.setField(selectionState || readSelectionState(form), "leftEnzyme", target.value),
        "Left enzyme updated locally. Click Preview Assembly to validate.",
        "text-warning"
      );
    } else if (target.name === "right_enzyme") {
      applySelectionUpdate(
        form,
        selectionStore.setField(selectionState || readSelectionState(form), "rightEnzyme", target.value),
        "Right enzyme updated locally. Click Preview Assembly to validate.",
        "text-warning"
      );
    } else if (target.name === "selected_enzymes") {
      clearFragmentAnchors();
      applySelectionUpdate(
        form,
        selectionStore.setSelectedEnzymes(selectionState || readSelectionState(form), selectedValues(target)),
        "Restriction enzyme selection updated. Select fragments from the map, then click Preview Assembly.",
        "text-warning"
      );
    } else if (target.name === "vector_fragment_index") {
      applySelectionUpdate(
        form,
        selectionStore.setField(selectionState || readSelectionState(form), "vectorFragmentIndex", target.value),
        "Vector fragment updated locally. Click Preview Assembly to validate.",
        "text-warning"
      );
    } else if (target.name === "insert_fragment_index") {
      applySelectionUpdate(
        form,
        selectionStore.setField(selectionState || readSelectionState(form), "insertFragmentIndex", target.value),
        "Insert fragment updated locally. Click Preview Assembly to validate.",
        "text-warning"
      );
    } else if (target.name === "is_circular") {
      applySelectionUpdate(
        form,
        selectionStore.setField(selectionState || readSelectionState(form), "isCircular", target.value),
        "Topology updated locally. Click Preview Assembly to validate.",
        "text-warning"
      );
    }
  });

  document.addEventListener("input", function (event) {
    const search = event.target.closest("[data-cloning-enzyme-search]");
    if (!search) {
      return;
    }
    enzymeTablePage = 1;
    refreshUnifiedEnzymeTable(getCloningPreviewForm());
  });

  document.addEventListener("click", function (event) {
    const previous = event.target.closest("[data-cloning-enzyme-prev-page]");
    const next = event.target.closest("[data-cloning-enzyme-next-page]");
    if (!previous && !next) {
      return;
    }
    if (previous) {
      enzymeTablePage = Math.max(1, enzymeTablePage - 1);
    } else {
      enzymeTablePage += 1;
    }
    refreshUnifiedEnzymeTable(getCloningPreviewForm());
  });

  document.addEventListener("change", function (event) {
    const checkbox = event.target.closest("[data-cloning-enzyme-toggle]");
    if (!checkbox) {
      return;
    }
    const form = getCloningPreviewForm();
    const row = checkbox.closest("[data-cloning-enzyme-row]");
    const enzymeName = row ? row.dataset.enzymeName || "" : "";
    if (!form || !enzymeName) {
      return;
    }
    const current = selectedValues(form.elements.selected_enzymes);
    const next = checkbox.checked
      ? current.concat([enzymeName]).filter((value, index, list) => list.indexOf(value) === index)
      : current.filter((value) => value !== enzymeName);
    setSelectedEnzymeValues(form, next);
  });

  document.addEventListener("keydown", function (event) {
    if (event.key !== "Enter" && event.key !== " ") {
      return;
    }
    const trigger = event.target.closest("[data-cloning-enzyme-picker], [data-cloning-region-picker]");
    if (!trigger) {
      return;
    }
    event.preventDefault();
    trigger.dispatchEvent(new MouseEvent("click", { bubbles: true, shiftKey: event.shiftKey }));
  });

  document.addEventListener("submit", function (event) {
    if (event.target && event.target.id === "cloning-assembly-preview-form") {
      setFormLoading(event.target, "Updating assembly preview...");
    }
  });

  const initialForm = getCloningPreviewForm();
  selectionState = readSelectionState(initialForm);
  if (selectionStore && selectionState) {
    authoritativeSelectionKey = selectionStore.selectionKey(selectionState);
  }
  renderLocalSelection(initialForm);
})();
