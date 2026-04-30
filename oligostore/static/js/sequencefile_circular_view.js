(function () {
  const restriction = window.SequenceRestriction;

  function getElement(id) {
    return document.getElementById(id);
  }

  function escapeHtml(value) {
    return String(value)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;")
      .replace(/'/g, "&#39;");
  }

  const normalize = restriction.normalize;

  function polarToCartesian(cx, cy, radius, angleInDegrees) {
    const angleInRadians = ((angleInDegrees - 90) * Math.PI) / 180;
    return {
      x: cx + radius * Math.cos(angleInRadians),
      y: cy + radius * Math.sin(angleInRadians),
    };
  }

  function describeArc(cx, cy, radius, startAngle, endAngle) {
    const normalizedEnd = endAngle <= startAngle ? endAngle + 360 : endAngle;
    const span = Math.max(normalizedEnd - startAngle, 0.8);
    const actualEnd = startAngle + span;
    const start = polarToCartesian(cx, cy, radius, actualEnd);
    const end = polarToCartesian(cx, cy, radius, startAngle);
    const largeArcFlag = span > 180 ? "1" : "0";
    return `M ${start.x} ${start.y} A ${radius} ${radius} 0 ${largeArcFlag} 0 ${end.x} ${end.y}`;
  }

  function featureColor(feature) {
    if (isCdsFeature(feature)) {
      return "#0ea5e9";
    }
    if (isPrimerBindingFeature(feature)) {
      return "#2563eb";
    }
    return "#10b981";
  }

  function clampPosition(value, length) {
    const position = Number(value);
    if (!Number.isFinite(position)) {
      return 1;
    }
    return Math.min(Math.max(1, position), length);
  }

  function isPrimerBindingFeature(feature) {
    const typeValue = normalize(feature.type);
    const labelValue = normalize(feature.label);
    const noteValue = normalize(feature.note);
    return (
      typeValue.includes("primer") ||
      typeValue.includes("primer_bind") ||
      labelValue.includes("primer") ||
      noteValue.includes("primer") ||
      labelValue.includes("_fw") ||
      labelValue.includes("_rv") ||
      labelValue.includes("_rev")
    );
  }

  function isCdsFeature(feature) {
    return normalize(feature.type) === "cds";
  }

  function isPrimerOrMiscFeature(feature) {
    const typeValue = normalize(feature.type);
    return isPrimerBindingFeature(feature) || typeValue === "misc_feature" || typeValue === "misc_features" || typeValue === "misc features" || typeValue === "misc";
  }

  function isMiscFeature(feature) {
    const typeValue = normalize(feature.type);
    return typeValue === "misc_feature" || typeValue === "misc_features" || typeValue === "misc features" || typeValue === "misc";
  }

  function shouldDisplayFeature(app, feature) {
    if (isCdsFeature(feature)) {
      return app.state.showCdsFeatures;
    }
    if (isPrimerBindingFeature(feature)) {
      return app.state.showPrimerFeatures;
    }
    if (isMiscFeature(feature)) {
      return app.state.showMiscFeatures;
    }
    return false;
  }

  function updateFeatureToggleButtons(app) {
    [
      [app.els.toggleCdsFeaturesBtn, app.state.showCdsFeatures],
      [app.els.togglePrimerFeaturesBtn, app.state.showPrimerFeatures],
      [app.els.toggleMiscFeaturesBtn, app.state.showMiscFeatures],
    ].forEach(([button, enabled]) => {
      if (!button) {
        return;
      }
      button.classList.toggle("btn-primary", enabled);
      button.classList.toggle("btn-outline", !enabled);
      button.setAttribute("aria-pressed", enabled ? "true" : "false");
    });
  }

  function circularFeatureRadius(feature, baseRadius) {
    return isCdsFeature(feature) ? baseRadius : baseRadius - 22;
  }

  function createApp(raw, config) {
    const recordSummaries = JSON.parse(raw.textContent);
    if (!recordSummaries.length) {
      return null;
    }

    return {
      recordSummaries,
      recordDetails: new Map(),
      urls: {
        recordDataUrl: config.dataset.recordDataUrl,
      },
      els: {
        recordSelect: getElement("circularRecordSelect"),
        zoomSlider: getElement("circularZoomSlider"),
        zoomInBtn: getElement("circularZoomInBtn"),
        zoomOutBtn: getElement("circularZoomOutBtn"),
        resetBtn: getElement("circularResetBtn"),
        toggleCdsFeaturesBtn: getElement("circularToggleCdsFeaturesBtn"),
        togglePrimerFeaturesBtn: getElement("circularTogglePrimerFeaturesBtn"),
        toggleMiscFeaturesBtn: getElement("circularToggleMiscFeaturesBtn"),
        zoomHint: getElement("circularZoomHint"),
        windowLabel: getElement("circularWindowLabel"),
        featureArcs: getElement("circularFeatureArcs"),
        restrictionMarks: getElement("circularRestrictionMarks"),
        axisTicks: getElement("circularAxisTicks"),
        centerTitle: getElement("circularCenterTitle"),
        centerMeta: getElement("circularCenterMeta"),
        centerSecondary: getElement("circularCenterSecondary"),
        mapStatus: getElement("circularMapStatus"),
        recordLength: getElement("circularRecordLength"),
        featureTotal: getElement("circularFeatureTotal"),
        restrictionCount: getElement("circularRestrictionCount"),
        clearSelectionBtn: getElement("circularClearSelectionBtn"),
        featureDetails: getElement("circularFeatureDetails"),
        selectionSummary: getElement("circularSelectionSummary"),
        featureSearchInput: getElement("circularFeatureSearchInput"),
        featureCount: getElement("circularFeatureCount"),
        featureTableBody: getElement("circularFeatureTableBody"),
        restrictionSearchInput: getElement("circularRestrictionSearchInput"),
        selectTopRestrictionEnzymesBtn: getElement("circularSelectTopRestrictionEnzymesBtn"),
        clearRestrictionSelectionBtn: getElement("circularClearRestrictionSelectionBtn"),
        restrictionSelectionCount: getElement("circularRestrictionSelectionCount"),
        restrictionTableBody: getElement("circularRestrictionTableBody"),
        restrictionPageInfo: getElement("circularRestrictionPageInfo"),
        restrictionPrevPageBtn: getElement("circularRestrictionPrevPageBtn"),
        restrictionNextPageBtn: getElement("circularRestrictionNextPageBtn"),
        restrictionSummary: getElement("circularRestrictionSummary"),
      },
      state: {
        recordIndex: 0,
        selectedFeatureIndex: null,
        featureQuery: "",
        loadError: "",
        loading: false,
        zoomPercent: 100,
        selectedRestrictionEnzymes: new Set(),
        restrictionEnzymeQuery: "",
        restrictionTablePage: 1,
        restrictionTablePageSize: 10,
        showCdsFeatures: true,
        showPrimerFeatures: true,
        showMiscFeatures: true,
      },
    };
  }

  function getRecord(app) {
    return app.recordDetails.get(app.state.recordIndex) || null;
  }

  async function ensureRecordLoaded(app) {
    const summary = app.recordSummaries[app.state.recordIndex];
    if (!summary || app.recordDetails.has(app.state.recordIndex)) {
      return;
    }

    app.state.loading = true;
    app.state.loadError = "";
    render(app);
    try {
      const response = await fetch(
        `${app.urls.recordDataUrl}?record_index=${app.state.recordIndex}&start=1&end=${summary.length}`
      );
      if (!response.ok) {
        throw new Error(`HTTP ${response.status}`);
      }
      app.recordDetails.set(app.state.recordIndex, await response.json());
    } catch (error) {
      app.state.loadError = `Could not load record data (${error.message}).`;
    } finally {
      app.state.loading = false;
      render(app);
    }
  }

  function getFilteredFeatures(app, record) {
    const query = normalize(app.state.featureQuery);
    if (!query) {
      return record.features.map((feature, index) => ({ feature, index }));
    }
    return record.features
      .map((feature, index) => ({ feature, index }))
      .filter(({ feature }) => {
        const haystack = [
          feature.label,
          feature.type,
          feature.source || "imported",
          feature.start,
          feature.end,
          feature.strand,
        ].map(normalize).join(" ");
        return haystack.includes(query);
      });
  }

  function getSortedRestrictionEnzymes(record) {
    return restriction.getSortedEnzymes(record);
  }

  function applyTopRestrictionSelection(app, record, limit = 5) {
    app.state.selectedRestrictionEnzymes = new Set();
    restriction.applyTopSelection(record, app.state.selectedRestrictionEnzymes, limit);
  }

  function renderFeatureDetails(app, feature) {
    if (!feature) {
      app.els.featureDetails.innerHTML = "No feature selected.";
      app.els.selectionSummary.textContent = "No feature selected on map.";
      return;
    }

    const sourceLabel = feature.source === "user" ? "User-added" : "Imported";
    const strandLabel = feature.strand === -1 ? "Reverse" : feature.strand === 1 ? "Forward" : "Unknown";
    const length = Number(feature.end) - Number(feature.start) + 1;
    app.els.featureDetails.innerHTML = `
      <div class="space-y-2">
        <div class="text-base font-semibold">${escapeHtml(feature.label || feature.type || "Feature")}</div>
        <div><span class="opacity-60">Type:</span> ${escapeHtml(feature.type || "-")}</div>
        <div><span class="opacity-60">Range:</span> ${Number(feature.start).toLocaleString()}-${Number(feature.end).toLocaleString()} (${length.toLocaleString()} bp)</div>
        <div><span class="opacity-60">Strand:</span> ${strandLabel}</div>
        <div><span class="opacity-60">Source:</span> ${sourceLabel}</div>
        <div class="opacity-80">${escapeHtml(feature.description || feature.note || "No description.")}</div>
      </div>
    `;
    app.els.selectionSummary.textContent = `${escapeHtml(feature.label || feature.type || "Feature")} selected at ${Number(feature.start).toLocaleString()}-${Number(feature.end).toLocaleString()}.`;
  }

  function renderFeatureTable(app, record) {
    const filtered = getFilteredFeatures(app, record);
    app.els.featureCount.textContent = `${filtered.length.toLocaleString()} feature(s)`;
    if (!filtered.length) {
      app.els.featureTableBody.innerHTML = '<tr><td colspan="6" class="text-center opacity-70 py-3">No features match the current filter.</td></tr>';
      return;
    }

    app.els.featureTableBody.innerHTML = filtered
      .map(({ feature, index }) => {
        const isSelected = app.state.selectedFeatureIndex === index;
        const sourceLabel = feature.source === "user" ? "User" : "Imported";
        return `
          <tr class="circular-feature-row${isSelected ? " is-selected" : ""}" data-feature-index="${index}">
            <td>${escapeHtml(feature.label || "-")}</td>
            <td>${escapeHtml(feature.type || "-")}</td>
            <td>${Number(feature.start).toLocaleString()}</td>
            <td>${Number(feature.end).toLocaleString()}</td>
            <td>${(Number(feature.end) - Number(feature.start) + 1).toLocaleString()}</td>
            <td>${sourceLabel}</td>
          </tr>
        `;
      })
      .join("");

    app.els.featureTableBody.querySelectorAll("[data-feature-index]").forEach((row) => {
      row.addEventListener("click", () => {
        app.state.selectedFeatureIndex = Number(row.dataset.featureIndex);
        render(app);
      });
    });
  }

  function renderRestrictionTable(app, record) {
    const sorted = getSortedRestrictionEnzymes(record);
    const available = new Set(sorted.map((item) => item.enzyme));
    app.state.selectedRestrictionEnzymes = new Set(
      [...app.state.selectedRestrictionEnzymes].filter((enzyme) => available.has(enzyme))
    );

    const query = normalize(app.state.restrictionEnzymeQuery);
    const selectedItems = sorted.filter((item) => app.state.selectedRestrictionEnzymes.has(item.enzyme));
    const unselectedItems = sorted.filter((item) => !app.state.selectedRestrictionEnzymes.has(item.enzyme));
    const filteredUnselectedItems = !query
      ? unselectedItems
      : unselectedItems.filter(({ enzyme }) => normalize(enzyme).includes(query));

    if (!sorted.length) {
      app.els.restrictionTableBody.innerHTML = '<tr><td colspan="5" class="text-center opacity-70 py-3">No restriction sites in this record.</td></tr>';
      app.els.restrictionSelectionCount.textContent = "0 selected";
      app.els.restrictionPageInfo.textContent = "Page 0 of 0";
      app.els.restrictionPrevPageBtn.disabled = true;
      app.els.restrictionNextPageBtn.disabled = true;
      app.els.restrictionSummary.textContent = "No restriction sites are available for this record.";
      return;
    }

    if (!selectedItems.length && !filteredUnselectedItems.length) {
      app.els.restrictionTableBody.innerHTML = '<tr><td colspan="5" class="text-center opacity-70 py-3">No enzymes match this search.</td></tr>';
      app.els.restrictionSelectionCount.textContent = `${app.state.selectedRestrictionEnzymes.size.toLocaleString()} selected | 0 matching`;
      app.els.restrictionPageInfo.textContent = "Page 0 of 0";
      app.els.restrictionPrevPageBtn.disabled = true;
      app.els.restrictionNextPageBtn.disabled = true;
      app.els.restrictionSummary.textContent = "No enzymes match the current filter.";
      return;
    }

    const totalPages = Math.max(1, Math.ceil(filteredUnselectedItems.length / app.state.restrictionTablePageSize));
    app.state.restrictionTablePage = Math.min(Math.max(1, app.state.restrictionTablePage), totalPages);
    const pageStart = (app.state.restrictionTablePage - 1) * app.state.restrictionTablePageSize;
    const pageItems = filteredUnselectedItems.slice(pageStart, pageStart + app.state.restrictionTablePageSize);
    const rowsToRender = [...selectedItems, ...pageItems];

    app.els.restrictionTableBody.innerHTML = "";
    rowsToRender.forEach(({ enzyme, count, site, cutOffset }) => {
      const row = document.createElement("tr");
      row.className = "hover";

      const input = document.createElement("input");
      input.type = "checkbox";
      input.className = "checkbox checkbox-xs checkbox-primary";
      input.value = enzyme;
      input.checked = app.state.selectedRestrictionEnzymes.has(enzyme);
      input.addEventListener("change", () => {
        if (input.checked) {
          app.state.selectedRestrictionEnzymes.add(enzyme);
        } else {
          app.state.selectedRestrictionEnzymes.delete(enzyme);
        }
        render(app);
      });

      const selectCell = document.createElement("td");
      selectCell.appendChild(input);

      const markerCell = document.createElement("td");
      markerCell.appendChild(restriction.createColorSwatch(enzyme));

      const enzymeCell = document.createElement("td");
      enzymeCell.className = "font-mono";
      enzymeCell.textContent = enzyme;

      const siteCell = document.createElement("td");
      siteCell.className = "font-mono";
      siteCell.textContent = restriction.formatCutSite(site, cutOffset);
      siteCell.title = `${site || "-"} | cut offset: ${Number.isFinite(cutOffset) ? cutOffset : "n/a"}`;

      const countCell = document.createElement("td");
      countCell.textContent = count.toLocaleString();

      row.appendChild(selectCell);
      row.appendChild(markerCell);
      row.appendChild(enzymeCell);
      row.appendChild(siteCell);
      row.appendChild(countCell);
      app.els.restrictionTableBody.appendChild(row);
    });

    const matchingTotal = selectedItems.length + filteredUnselectedItems.length;
    app.els.restrictionSelectionCount.textContent = `${app.state.selectedRestrictionEnzymes.size.toLocaleString()} selected | ${matchingTotal.toLocaleString()} matching`;
    app.els.restrictionPageInfo.textContent = `Page ${app.state.restrictionTablePage.toLocaleString()} of ${totalPages.toLocaleString()}`;
    app.els.restrictionPrevPageBtn.disabled = app.state.restrictionTablePage <= 1;
    app.els.restrictionNextPageBtn.disabled = app.state.restrictionTablePage >= totalPages;
    app.els.restrictionSummary.textContent = app.state.selectedRestrictionEnzymes.size
      ? `${app.state.selectedRestrictionEnzymes.size.toLocaleString()} selected enzyme(s) highlighted on the circular map.`
      : "Select one or more enzymes from the table to highlight restriction sites on the map.";
  }

  function renderCircularMap(app, record) {
    const cx = 350;
    const cy = 350;
    const zoomScale = app.state.zoomPercent / 100;
    const zoomOffset = (zoomScale - 1) * 42;
    const axisRadius = 255 + zoomOffset;
    const featureBaseRadius = 230 + zoomOffset;
    const restrictionInner = 268 + zoomOffset;
    const restrictionOuter = 288 + zoomOffset;
    const length = Math.max(1, Number(record.length) || 1);
    const selectedFeature = app.state.selectedFeatureIndex === null ? null : record.features[app.state.selectedFeatureIndex];
    const selectedFeatureVisible = selectedFeature ? shouldDisplayFeature(app, selectedFeature) : false;
    const displayFeatures = record.features
      .map((feature, index) => ({ feature, index }))
      .filter(({ feature }) => shouldDisplayFeature(app, feature));

    app.els.axisTicks.innerHTML = "";
    app.els.featureArcs.innerHTML = "";
    app.els.restrictionMarks.innerHTML = "";

    const tickCount = Math.min(24, Math.max(8, Math.ceil(length / 500)));
    for (let index = 0; index < tickCount; index += 1) {
      const position = 1 + Math.round((index / tickCount) * (length - 1));
      const angle = (position / length) * 360;
      const inner = polarToCartesian(cx, cy, axisRadius - 10, angle);
      const outer = polarToCartesian(cx, cy, axisRadius + 10, angle);
      const textPoint = polarToCartesian(cx, cy, axisRadius + 28, angle);

      const line = document.createElementNS("http://www.w3.org/2000/svg", "line");
      line.setAttribute("x1", inner.x);
      line.setAttribute("y1", inner.y);
      line.setAttribute("x2", outer.x);
      line.setAttribute("y2", outer.y);
      line.setAttribute("stroke", "#94a3b8");
      line.setAttribute("stroke-width", "1.5");
      app.els.axisTicks.appendChild(line);

      const label = document.createElementNS("http://www.w3.org/2000/svg", "text");
      label.setAttribute("x", textPoint.x);
      label.setAttribute("y", textPoint.y);
      label.setAttribute("text-anchor", "middle");
      label.setAttribute("fill", "#64748b");
      label.setAttribute("font-size", "11");
      label.textContent = position.toLocaleString();
      app.els.axisTicks.appendChild(label);
    }

    displayFeatures.forEach(({ feature, index }) => {
      const start = clampPosition(feature.start, length);
      const end = clampPosition(feature.end, length);
      const startAngle = ((start - 1) / length) * 360;
      const endAngle = (end / length) * 360;
      const radius = circularFeatureRadius(feature, featureBaseRadius);
      const isSelected = selectedFeature === feature;
      const path = document.createElementNS("http://www.w3.org/2000/svg", "path");
      path.setAttribute("d", describeArc(cx, cy, radius, startAngle, endAngle));
      path.setAttribute("fill", "none");
      path.setAttribute("stroke", featureColor(feature));
      path.setAttribute("stroke-width", isSelected ? "13" : "8");
      path.setAttribute("stroke-linecap", "round");
      path.setAttribute("opacity", selectedFeatureVisible && !isSelected ? "0.32" : "0.96");
      path.style.cursor = "pointer";
      path.addEventListener("click", () => {
        app.state.selectedFeatureIndex = index;
        render(app);
      });

      const title = document.createElementNS("http://www.w3.org/2000/svg", "title");
      title.textContent = `${feature.label || feature.type || "Feature"} (${feature.start}-${feature.end})`;
      path.appendChild(title);
      app.els.featureArcs.appendChild(path);

      const spanDegrees = Math.abs(endAngle - startAngle);
      if (isSelected || spanDegrees >= 24) {
        const labelAngle = startAngle + (spanDegrees / 2);
        const labelPoint = polarToCartesian(cx, cy, radius - 18, labelAngle);
        const label = document.createElementNS("http://www.w3.org/2000/svg", "text");
        label.setAttribute("x", labelPoint.x);
        label.setAttribute("y", labelPoint.y);
        label.setAttribute("text-anchor", "middle");
        label.setAttribute("dominant-baseline", "middle");
        label.setAttribute("fill", isSelected ? "#0f172a" : "#334155");
        label.setAttribute("font-size", isSelected ? "12" : "10");
        label.setAttribute("font-weight", isSelected ? "700" : "600");
        label.setAttribute("paint-order", "stroke");
        label.setAttribute("stroke", "#ffffff");
        label.setAttribute("stroke-width", "4");
        label.setAttribute("pointer-events", "none");
        label.textContent = String(feature.label || feature.type || "Feature").slice(0, 22);
        app.els.featureArcs.appendChild(label);
      }
    });

    const hasSelection = app.state.selectedRestrictionEnzymes.size > 0;
    const sitesToRender = hasSelection
      ? record.restriction_sites.filter((site) => app.state.selectedRestrictionEnzymes.has(site.enzyme))
      : [];
    const selectedEnzymes = [...app.state.selectedRestrictionEnzymes];
    sitesToRender.forEach((site) => {
      const angle = ((Number(site.start) - 1) / length) * 360;
      const enzymeLane = Math.max(0, selectedEnzymes.indexOf(site.enzyme));
      const laneOffset = Math.min(enzymeLane, 4) * 5;
      const inner = polarToCartesian(cx, cy, restrictionInner + laneOffset, angle);
      const outer = polarToCartesian(cx, cy, restrictionOuter + laneOffset, angle);
      const color = restriction.getEnzymeColor(site.enzyme);
      const mark = document.createElementNS("http://www.w3.org/2000/svg", "line");
      mark.setAttribute("x1", inner.x);
      mark.setAttribute("y1", inner.y);
      mark.setAttribute("x2", outer.x);
      mark.setAttribute("y2", outer.y);
      mark.setAttribute("stroke", color);
      mark.setAttribute("stroke-width", hasSelection ? "3.5" : "2.5");
      mark.setAttribute("stroke-linecap", "round");
      mark.setAttribute("opacity", hasSelection ? "0.95" : "0.55");

      const title = document.createElementNS("http://www.w3.org/2000/svg", "title");
      title.textContent = restriction.siteTitle(site);
      mark.appendChild(title);
      app.els.restrictionMarks.appendChild(mark);

      if (sitesToRender.length <= 40) {
        const labelPoint = polarToCartesian(cx, cy, restrictionOuter + laneOffset + 18, angle);
        const label = document.createElementNS("http://www.w3.org/2000/svg", "text");
        label.setAttribute("x", labelPoint.x);
        label.setAttribute("y", labelPoint.y);
        label.setAttribute("text-anchor", "middle");
        label.setAttribute("dominant-baseline", "middle");
        label.setAttribute("fill", color);
        label.setAttribute("font-size", "10");
        label.setAttribute("font-weight", "700");
        label.setAttribute("paint-order", "stroke");
        label.setAttribute("stroke", "#ffffff");
        label.setAttribute("stroke-width", "3");
        label.setAttribute("pointer-events", "none");
        label.textContent = String(site.enzyme || "").slice(0, 8);
        app.els.restrictionMarks.appendChild(label);
      }
    });
  }

  function render(app) {
    const summary = app.recordSummaries[app.state.recordIndex];
    const record = getRecord(app);
    updateFeatureToggleButtons(app);
    app.els.centerTitle.textContent = summary ? summary.id : "Record";
    app.els.centerMeta.textContent = summary ? `${Number(summary.length).toLocaleString()} bp` : "-";
    app.els.centerSecondary.textContent = "Circular sequence map";
    app.els.recordLength.textContent = summary ? `${Number(summary.length).toLocaleString()} bp` : "-";
    app.els.zoomSlider.value = String(app.state.zoomPercent);
    app.els.zoomHint.textContent = `${app.state.zoomPercent}% ring scale`;
    app.els.windowLabel.textContent = summary ? `${summary.id} | full record` : "Full record";

    if (!record) {
      app.els.mapStatus.textContent = app.state.loading ? "Loading circular map..." : (app.state.loadError || "Record data not loaded.");
      app.els.featureTotal.textContent = "-";
      app.els.restrictionCount.textContent = "-";
      app.els.featureTableBody.innerHTML = '<tr><td colspan="6" class="text-center opacity-70 py-3">No feature data loaded.</td></tr>';
      app.els.featureCount.textContent = "";
      app.els.restrictionTableBody.innerHTML = '<tr><td colspan="5" class="text-center opacity-70 py-3">No restriction data loaded.</td></tr>';
      app.els.restrictionSelectionCount.textContent = "0 selected";
      app.els.restrictionPageInfo.textContent = "Page 0 of 0";
      app.els.restrictionSummary.textContent = app.state.loading ? "Loading restriction sites..." : "Restriction data not loaded.";
      app.els.axisTicks.innerHTML = "";
      app.els.featureArcs.innerHTML = "";
      app.els.restrictionMarks.innerHTML = "";
      renderFeatureDetails(app, null);
      return;
    }

    app.els.featureTotal.textContent = `${record.features.length.toLocaleString()} total`;
    app.els.restrictionCount.textContent = `${record.restriction_sites.length.toLocaleString()} total`;
    const displayedFeatureCount = record.features.filter((feature) => shouldDisplayFeature(app, feature)).length;
    app.els.mapStatus.textContent = `${displayedFeatureCount.toLocaleString()} displayed feature(s): CDS and primer_bind/misc_feature lanes | ${record.restriction_sites.length.toLocaleString()} restriction site(s) available`;
    renderCircularMap(app, record);
    renderFeatureTable(app, record);
    renderRestrictionTable(app, record);
    renderFeatureDetails(app, app.state.selectedFeatureIndex === null ? null : record.features[app.state.selectedFeatureIndex]);
  }

  document.addEventListener("DOMContentLoaded", () => {
    const raw = getElement("sequence-circular-records-data");
    const config = getElement("sequence-circular-config");
    if (!raw || !config) {
      return;
    }

    const app = createApp(raw, config);
    if (!app) {
      return;
    }

    app.recordSummaries.forEach((record, index) => {
      const option = document.createElement("option");
      option.value = String(index);
      option.textContent = `${record.id} (${Number(record.length).toLocaleString()} bp)`;
      app.els.recordSelect.appendChild(option);
    });

    app.els.recordSelect.addEventListener("change", () => {
      app.state.recordIndex = Number(app.els.recordSelect.value);
      app.state.selectedFeatureIndex = null;
      app.state.featureQuery = "";
      app.state.selectedRestrictionEnzymes = new Set();
      app.state.restrictionEnzymeQuery = "";
      app.state.restrictionTablePage = 1;
      app.els.featureSearchInput.value = "";
      app.els.restrictionSearchInput.value = "";
      render(app);
      void ensureRecordLoaded(app);
    });
    app.els.zoomSlider.addEventListener("input", () => {
      app.state.zoomPercent = Number(app.els.zoomSlider.value);
      render(app);
    });
    app.els.zoomInBtn.addEventListener("click", () => {
      app.state.zoomPercent = Math.min(130, app.state.zoomPercent + 10);
      render(app);
    });
    app.els.zoomOutBtn.addEventListener("click", () => {
      app.state.zoomPercent = Math.max(70, app.state.zoomPercent - 10);
      render(app);
    });
    app.els.resetBtn.addEventListener("click", () => {
      app.state.zoomPercent = 100;
      render(app);
    });
    app.els.toggleCdsFeaturesBtn.addEventListener("click", () => {
      app.state.showCdsFeatures = !app.state.showCdsFeatures;
      render(app);
    });
    app.els.togglePrimerFeaturesBtn.addEventListener("click", () => {
      app.state.showPrimerFeatures = !app.state.showPrimerFeatures;
      render(app);
    });
    app.els.toggleMiscFeaturesBtn.addEventListener("click", () => {
      app.state.showMiscFeatures = !app.state.showMiscFeatures;
      render(app);
    });
    app.els.clearSelectionBtn.addEventListener("click", () => {
      app.state.selectedFeatureIndex = null;
      render(app);
    });
    app.els.featureSearchInput.addEventListener("input", () => {
      app.state.featureQuery = app.els.featureSearchInput.value;
      render(app);
    });
    app.els.restrictionSearchInput.addEventListener("input", () => {
      app.state.restrictionEnzymeQuery = app.els.restrictionSearchInput.value;
      app.state.restrictionTablePage = 1;
      render(app);
    });
    app.els.selectTopRestrictionEnzymesBtn.addEventListener("click", () => {
      const record = getRecord(app);
      if (!record) {
        return;
      }
      applyTopRestrictionSelection(app, record);
      app.state.restrictionTablePage = 1;
      render(app);
    });
    app.els.clearRestrictionSelectionBtn.addEventListener("click", () => {
      app.state.selectedRestrictionEnzymes = new Set();
      app.state.restrictionTablePage = 1;
      render(app);
    });
    app.els.restrictionPrevPageBtn.addEventListener("click", () => {
      app.state.restrictionTablePage = Math.max(1, app.state.restrictionTablePage - 1);
      render(app);
    });
    app.els.restrictionNextPageBtn.addEventListener("click", () => {
      app.state.restrictionTablePage += 1;
      render(app);
    });

    app.els.recordSelect.value = "0";
    render(app);
    void ensureRecordLoaded(app);
  });
})();
