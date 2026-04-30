(function () {
  const ns = window.SequenceLinearView = window.SequenceLinearView || {};
  const restriction = window.SequenceRestriction;

  ns.registerFeatureModule = function registerFeatureModule(app) {
    const { state, els } = app;

    app.normalize = function normalize(value) {
      return String(value || "").toLowerCase();
    };

    app.getRestrictionHitsInWindow = function getRestrictionHitsInWindow(record, start, end) {
      return (record.restriction_sites || []).filter((site) => !(site.end < start || site.start > end));
    };

    app.getSortedRestrictionEnzymes = function getSortedRestrictionEnzymes(record) {
      return restriction.getSortedEnzymes(record);
    };

    app.applyTopRestrictionSelection = function applyTopRestrictionSelection(record, limit = 5) {
      state.selectedRestrictionEnzymes = new Set();
      restriction.applyTopSelection(record, state.selectedRestrictionEnzymes, limit);
    };

    app.renderRestrictionEnzymeTable = function renderRestrictionEnzymeTable(record) {
      const sorted = app.getSortedRestrictionEnzymes(record);
      const available = new Set(sorted.map((item) => item.enzyme));
      state.selectedRestrictionEnzymes = new Set(
        [...state.selectedRestrictionEnzymes].filter((enzyme) => available.has(enzyme))
      );
      const query = app.normalize(state.restrictionEnzymeQuery).trim();
      const selectedItems = sorted.filter((item) => state.selectedRestrictionEnzymes.has(item.enzyme));
      const unselectedItems = sorted.filter((item) => !state.selectedRestrictionEnzymes.has(item.enzyme));
      const filteredUnselectedItems = !query
        ? unselectedItems
        : unselectedItems.filter(({ enzyme }) => app.normalize(enzyme).includes(query));

      if (!sorted.length) {
        els.restrictionEnzymeTableBody.innerHTML = '<tr><td colspan="5" class="text-center opacity-70 py-3">No restriction sites in this record.</td></tr>';
        els.restrictionSelectionCount.textContent = "0 selected";
        els.restrictionEnzymePageInfo.textContent = "Page 0 of 0";
        els.restrictionEnzymePrevPageBtn.disabled = true;
        els.restrictionEnzymeNextPageBtn.disabled = true;
        return;
      }

      if (!selectedItems.length && !filteredUnselectedItems.length) {
        els.restrictionEnzymeTableBody.innerHTML = '<tr><td colspan="5" class="text-center opacity-70 py-3">No enzymes match this search.</td></tr>';
        els.restrictionSelectionCount.textContent = `${state.selectedRestrictionEnzymes.size.toLocaleString()} selected | 0 matching`;
        els.restrictionEnzymePageInfo.textContent = "Page 0 of 0";
        els.restrictionEnzymePrevPageBtn.disabled = true;
        els.restrictionEnzymeNextPageBtn.disabled = true;
        return;
      }

      const totalPages = Math.max(1, Math.ceil(filteredUnselectedItems.length / state.restrictionTablePageSize));
      state.restrictionTablePage = Math.min(Math.max(1, state.restrictionTablePage), totalPages);
      const pageStart = (state.restrictionTablePage - 1) * state.restrictionTablePageSize;
      const pageEnd = pageStart + state.restrictionTablePageSize;
      const pageItems = filteredUnselectedItems.slice(pageStart, pageEnd);
      const rowsToRender = [...selectedItems, ...pageItems];

      els.restrictionEnzymeTableBody.innerHTML = "";
      rowsToRender.forEach(({ enzyme, count, site, cutOffset }) => {
        const row = document.createElement("tr");
        row.className = "hover";
        const input = document.createElement("input");
        input.type = "checkbox";
        input.className = "checkbox checkbox-xs checkbox-primary";
        input.value = enzyme;
        input.checked = state.selectedRestrictionEnzymes.has(enzyme);
        input.addEventListener("change", () => {
          if (input.checked) {
            state.selectedRestrictionEnzymes.add(enzyme);
          } else {
            state.selectedRestrictionEnzymes.delete(enzyme);
          }
          app.render();
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
        els.restrictionEnzymeTableBody.appendChild(row);
      });

      const matchingTotal = selectedItems.length + filteredUnselectedItems.length;
      els.restrictionSelectionCount.textContent = `${state.selectedRestrictionEnzymes.size.toLocaleString()} selected | ${matchingTotal.toLocaleString()} matching`;
      els.restrictionEnzymePageInfo.textContent = `Page ${state.restrictionTablePage.toLocaleString()} of ${totalPages.toLocaleString()}`;
      els.restrictionEnzymePrevPageBtn.disabled = state.restrictionTablePage <= 1;
      els.restrictionEnzymeNextPageBtn.disabled = state.restrictionTablePage >= totalPages;
    };

    app.getFilteredFeatureIndexes = function getFilteredFeatureIndexes(record) {
      const query = app.normalize(state.featureQuery).trim();
      if (!query) {
        return record.features.map((_, index) => index);
      }
      return record.features
        .map((feature, index) => ({ feature, index }))
        .filter(({ feature }) => {
          const haystack = [
            feature.label,
            feature.type,
            feature.start,
            feature.end,
            feature.strand,
          ].map(app.normalize).join(" ");
          return haystack.includes(query);
        })
        .map(({ index }) => index);
    };

    app.getTableFeatureIndexes = function getTableFeatureIndexes(record) {
      const filtered = app.getFilteredFeatureIndexes(record);
      if (!state.showMapSelectedOnly) {
        return filtered;
      }
      return filtered.filter((index) => state.mapSelectedFeatureIndexes.has(index));
    };

    app.featureBounds = function featureBounds(feature) {
      const rawStart = Number(feature.start);
      const rawEnd = Number(feature.end);
      const safeStart = Number.isFinite(rawStart) ? rawStart : 1;
      const safeEnd = Number.isFinite(rawEnd) ? rawEnd : safeStart;
      return {
        start: Math.min(safeStart, safeEnd),
        end: Math.max(safeStart, safeEnd),
      };
    };

    app.featureLength = function featureLength(feature) {
      const bounds = app.featureBounds(feature);
      return bounds.end - bounds.start + 1;
    };

    app.featureDescription = function featureDescription(feature) {
      const description = String(
        feature.description || feature.note || feature.product || feature.label || feature.type || ""
      ).trim();
      return description || "No description available.";
    };

    app.isPrimerBindingFeature = function isPrimerBindingFeature(feature) {
      const typeValue = app.normalize(feature.type);
      const labelValue = app.normalize(feature.label);
      const noteValue = app.normalize(feature.note);
      if (typeValue.includes("primer") || typeValue.includes("primer_bind")) {
        return true;
      }
      if (labelValue.includes("primer") || noteValue.includes("primer")) {
        return true;
      }
      return labelValue.includes("_fw") || labelValue.includes("_rv") || labelValue.includes("_rev");
    };

    app.isCdsFeature = function isCdsFeature(feature) {
      return app.normalize(feature.type) === "cds";
    };

    app.isPrimerOrMiscFeature = function isPrimerOrMiscFeature(feature) {
      const typeValue = app.normalize(feature.type);
      return app.isPrimerBindingFeature(feature) || typeValue === "misc_feature" || typeValue === "misc_features" || typeValue === "misc features" || typeValue === "misc";
    };

    app.assignOverlapLane = function assignOverlapLane(laneEnds, visibleStart, visibleEnd, maxLanes = 5) {
      let lane = laneEnds.findIndex((end) => visibleStart > end);
      if (lane < 0) {
        lane = Math.min(laneEnds.length, maxLanes - 1);
      }
      laneEnds[lane] = Math.max(laneEnds[lane] || 0, visibleEnd);
      return lane;
    };

    app.normalizedStrandValue = function normalizedStrandValue(strand) {
      if (strand === -1 || strand === 1) {
        return strand;
      }
      const numeric = Number(strand);
      if (numeric === -1 || numeric === 1) {
        return numeric;
      }
      const strandText = app.normalize(strand);
      if (strandText.includes("reverse") || strandText.includes("complement") || strandText === "-") {
        return -1;
      }
      if (strandText.includes("forward") || strandText === "+") {
        return 1;
      }
      return 1;
    };

    app.strandLabel = function strandLabel(strand) {
      if (strand === 1) return "Forward (+)";
      if (strand === -1) return "Reverse (-)";
      return "Unknown";
    };

    app.jumpToFeature = function jumpToFeature(record, feature) {
      const length = app.featureLength(feature);
      const bounds = app.featureBounds(feature);
      const maxWindow = Math.min(5000, record.length);
      state.windowSize = Math.max(app.constants.MIN_WINDOW_BP, Math.min(maxWindow, length));
      const midpoint = Math.floor((bounds.start + bounds.end) / 2);
      const centeredStart = midpoint - Math.floor(state.windowSize / 2);
      const maxStart = Math.max(1, record.length - state.windowSize + 1);
      state.start = Math.max(1, Math.min(centeredStart, maxStart));
      const centeredMapStart = midpoint - Math.floor(state.mapWindowSize / 2);
      const maxMapStart = Math.max(1, record.length - state.mapWindowSize + 1);
      state.mapStart = Math.max(1, Math.min(centeredMapStart, maxMapStart));
    };

    app.applyFeatureSelection = function applyFeatureSelection(record, index, feature, isShift) {
      state.selectedFeatureIndex = index;
      if (isShift) {
        state.mapSelectedFeatureIndexes.add(index);
      } else {
        state.mapSelectedFeatureIndexes = new Set([index]);
      }
      state.showMapSelectedOnly = true;
      state.pendingFocusFeatureIndex = index;
      if (isShift && app.isPrimerBindingFeature(feature)) {
        if (state.selectedForwardPrimerIndex !== null && state.selectedReversePrimerIndex !== null) {
          state.selectedForwardPrimerIndex = null;
          state.selectedReversePrimerIndex = null;
          state.mapSelectedFeatureIndexes = new Set([index]);
        }
        const strand = app.normalizedStrandValue(feature.strand);
        if (strand === -1) {
          state.selectedReversePrimerIndex = index;
        } else {
          state.selectedForwardPrimerIndex = index;
        }
      }
      state.tablePage = 1;
      const didFocusPCRProduct = isShift && app.isPrimerBindingFeature(feature) && app.focusWindowOnPCRProduct(record);
      if (!didFocusPCRProduct) {
        app.jumpToFeature(record, feature);
      }
      app.render();
    };

    app.setWindowFeatureHoverText = function setWindowFeatureHoverText(feature) {
      if (!feature) {
        els.windowFeatureHoverDescription.textContent = "Hover a feature line to inspect its description.";
        return;
      }
      els.windowFeatureHoverDescription.textContent =
        `${feature.label} (${feature.type}) ${feature.start.toLocaleString()}-${feature.end.toLocaleString()} | ` +
        `${app.featureLength(feature).toLocaleString()} bp | ${app.strandLabel(feature.strand)} | ${app.featureDescription(feature)}`;
    };

    app.getHighlightedRestrictionPositions = function getHighlightedRestrictionPositions(windowRestrictionHits, start, end) {
      const highlightedPositions = new Map();
      windowRestrictionHits.forEach((site) => {
        const siteStart = Math.max(start, site.start);
        const siteEnd = Math.min(end, site.end);
        for (let position = siteStart; position <= siteEnd; position += 1) {
          if (!highlightedPositions.has(position)) {
            highlightedPositions.set(position, []);
          }
          highlightedPositions.get(position).push(site);
        }
      });
      return highlightedPositions;
    };

    app.renderStrandChunks = function renderStrandChunks(lineSegment, lineStart, highlightedPositions, strand, baseTransform = (base) => base) {
      const chunks = [];
      for (let offset = 0; offset < lineSegment.length; offset += app.constants.SEQUENCE_WINDOW_CHUNK_SIZE) {
        const chunk = lineSegment.slice(offset, offset + app.constants.SEQUENCE_WINDOW_CHUNK_SIZE);
        const chunkBaseIndex = lineStart + offset;
        const renderedChunk = chunk.split("").map((base, chunkOffset) => {
          const position = chunkBaseIndex + chunkOffset;
          return app.baseSpan(baseTransform(base), highlightedPositions.get(position) || [], position, strand);
        }).join("");
        chunks.push(renderedChunk);
      }
      return chunks.join(" ");
    };

    app.renderLineRestrictionSites = function renderLineRestrictionSites(windowRestrictionHits, lineStart, lineEnd) {
      return windowRestrictionHits
        .filter((site) => !(site.end < lineStart || site.start > lineEnd))
        .map((site) => `<span class="restriction-line-chip" style="--restriction-color: ${restriction.getEnzymeColor(site.enzyme)};">${app.escapeHtml(site.enzyme)} ${site.start.toLocaleString()}-${site.end.toLocaleString()}</span>`);
    };

    app.sequenceOffsetToDisplayOffset = function sequenceOffsetToDisplayOffset(offset) {
      return offset + Math.floor(offset / app.constants.SEQUENCE_WINDOW_CHUNK_SIZE);
    };

    app.sequenceDisplayLength = function sequenceDisplayLength(lineLength) {
      if (lineLength <= 0) {
        return 0;
      }
      return lineLength + Math.floor((lineLength - 1) / app.constants.SEQUENCE_WINDOW_CHUNK_SIZE);
    };

    app.buildRestrictionAnnotationRows = function buildRestrictionAnnotationRows(windowRestrictionHits, lineStart, lineEnd) {
      const lineLength = lineEnd - lineStart + 1;
      const displayLength = app.sequenceDisplayLength(lineLength);
      if (displayLength <= 0) {
        return null;
      }

      const markerChars = Array(displayLength).fill(" ");
      const labelChars = Array(displayLength).fill(" ");
      const lineSites = windowRestrictionHits.filter((site) => !(site.end < lineStart || site.start > lineEnd));

      lineSites.forEach((site) => {
        const startOffset = Math.max(0, site.start - lineStart);
        const cutOffset = Number.isFinite(Number(site.cut_offset)) ? Math.trunc(Number(site.cut_offset)) : 0;
        const cutPosition = site.start + cutOffset;
        const cutOffsetClamped = Math.max(lineStart, Math.min(lineEnd, cutPosition)) - lineStart;
        const markerIndex = app.sequenceOffsetToDisplayOffset(cutOffsetClamped);
        markerChars[markerIndex] = "v";

        const label = String(site.enzyme || "").trim();
        if (!label) {
          return;
        }
        const labelStart = app.sequenceOffsetToDisplayOffset(startOffset);
        for (let i = 0; i < label.length; i += 1) {
          const idx = labelStart + i;
          if (idx >= labelChars.length) {
            break;
          }
          if (labelChars[idx] === " ") {
            labelChars[idx] = label[i];
          }
        }
      });

      const markerText = markerChars.join("");
      const labelText = labelChars.join("");
      if (!markerText.trim() && !labelText.trim()) {
        return null;
      }
      return { markerText, labelText };
    };

    app.getLineNumberLayout = function getLineNumberLayout(end) {
      const width = String(Number(end).toLocaleString()).length;
      const prefix = " ".repeat(width + 2);
      return { width, prefix };
    };

    app.getSequenceSlice = function getSequenceSlice(record, start, end) {
      const regionStart = Number(record.region_start);
      const regionEnd = Number(record.region_end);
      if (!Number.isFinite(regionStart) || !Number.isFinite(regionEnd)) {
        return null;
      }
      if (start < regionStart || end > regionEnd) {
        return null;
      }
      const localStart = start - regionStart;
      const localEndExclusive = end - regionStart + 1;
      return record.sequence.slice(localStart, localEndExclusive);
    };

    app.renderSequenceWindow = function renderSequenceWindow(record, start, end) {
      const visible = app.getSequenceSlice(record, start, end);
      if (visible === null) {
        els.restrictionSiteSummary.textContent = state.loadingRecord
          ? "Loading sequence window..."
          : "Sequence window is outside the loaded region. Adjust pan/zoom to load it.";
        els.sequenceWindow.textContent = state.loadingRecord ? "Loading sequence..." : "Region not loaded.";
        return;
      }
      const lines = [];
      const hasRestrictionSelection = state.selectedRestrictionEnzymes.size > 0;
      const windowRestrictionHits = app.getRestrictionHitsInWindow(record, start, end)
        .filter((site) => state.selectedRestrictionEnzymes.has(site.enzyme));
      const highlightedPositions = hasRestrictionSelection
        ? app.getHighlightedRestrictionPositions(windowRestrictionHits, start, end)
        : new Map();
      const lineNumberLayout = app.getLineNumberLayout(end);

      for (let i = 0; i < visible.length; i += app.constants.SEQUENCE_WINDOW_LINE_SIZE) {
        const lineStart = start + i;
        const lineEnd = Math.min(end, lineStart + app.constants.SEQUENCE_WINDOW_LINE_SIZE - 1);
        const lineSegment = visible.slice(i, i + app.constants.SEQUENCE_WINDOW_LINE_SIZE);
        const forwardLine = app.renderStrandChunks(lineSegment, lineStart, highlightedPositions, "forward");
        const reverseLine = app.renderStrandChunks(lineSegment, lineStart, highlightedPositions, "reverse", app.complementBase);
        const lineNumber = lineStart.toLocaleString().padStart(lineNumberLayout.width, " ");

        lines.push(`${lineNumber}  5' ${forwardLine} 3'`);
        lines.push(`${lineNumberLayout.prefix}<span class="bottom-strand-row">3' ${reverseLine} 5'</span>`);
        if (hasRestrictionSelection) {
          const annotationRows = app.buildRestrictionAnnotationRows(windowRestrictionHits, lineStart, lineEnd);
          if (annotationRows) {
            lines.push(`${lineNumberLayout.prefix}   <span class="restriction-marker-row"><span class="restriction-cut-marker">${annotationRows.markerText}</span></span>`);
            lines.push(`${lineNumberLayout.prefix}   <span class="restriction-marker-row">${app.escapeHtml(annotationRows.labelText)}</span>`);
          }
          const lineSites = app.renderLineRestrictionSites(windowRestrictionHits, lineStart, lineEnd);
          if (lineSites.length) {
            lines.push(`${lineNumberLayout.prefix}<span class="restriction-row">RE: ${lineSites.join(" | ")}</span>`);
          }
        }
      }
      els.restrictionSiteSummary.textContent = hasRestrictionSelection
        ? `${windowRestrictionHits.length.toLocaleString()} restriction site(s) in current window`
        : "Select one or more enzymes from the table to visualize binding regions";
      els.sequenceWindow.innerHTML = lines.join("\n");
    };

    app.renderFeatureTrack = function renderFeatureTrack(record) {
      els.featureTrack.innerHTML = "";
      els.cdsFeatureTrack.innerHTML = "";
      els.primerMiscFeatureTrack.innerHTML = "";
      const visibleFeatureIndexes = new Set(app.getFilteredFeatureIndexes(record));
      const mapStart = state.mapStart;
      const mapEnd = Math.min(record.length, mapStart + state.mapWindowSize - 1);
      const mapLength = mapEnd - mapStart + 1;
      let cdsCount = 0;
      let primerMiscCount = 0;
      const cdsLaneEnds = [];
      const primerMiscLaneEnds = [];
      const laneHeight = 18;
      const laneGap = 4;

      function buildMarker(feature, index, color, shape) {
        const marker = document.createElement("div");
        marker.style.position = "absolute";
        marker.style.top = "0";
        marker.style.height = "100%";
        marker.style.backgroundColor = color;
        marker.style.border = "1px solid rgba(15, 23, 42, 0.35)";
        marker.style.borderRadius = "9999px";
        if (shape === "arrow-forward") {
          marker.style.clipPath = "polygon(0 0, calc(100% - 8px) 0, 100% 50%, calc(100% - 8px) 100%, 0 100%)";
        } else if (shape === "arrow-reverse") {
          marker.style.clipPath = "polygon(8px 0, 100% 0, 100% 100%, 8px 100%, 0 50%)";
        }
        marker.title = `${feature.label} (${feature.type}) ${feature.start}-${feature.end}`;
        marker.addEventListener("mousedown", (event) => handleFeatureMouseDown(event, index, feature));
        return marker;
      }

      function handleFeatureMouseDown(event, index, feature) {
        if (event.button !== 0) {
          return;
        }
        event.preventDefault();
        app.applyFeatureSelection(record, index, feature, event.shiftKey);
      }

      const visibleFeatures = record.features
        .map((feature, index) => ({ feature, index }))
        .filter(({ feature }) => {
          const bounds = app.featureBounds(feature);
          return !(bounds.end < mapStart || bounds.start > mapEnd);
        })
        .sort((a, b) => {
          const aBounds = app.featureBounds(a.feature);
          const bBounds = app.featureBounds(b.feature);
          const startDelta = aBounds.start - bBounds.start;
          if (startDelta !== 0) {
            return startDelta;
          }
          return (aBounds.end - aBounds.start) - (bBounds.end - bBounds.start);
        });

      visibleFeatures.forEach(({ feature, index }) => {
        const bounds = app.featureBounds(feature);
        const visibleStart = Math.max(bounds.start, mapStart);
        const visibleEnd = Math.min(bounds.end, mapEnd);
        const left = ((visibleStart - mapStart) / mapLength) * 100;
        const width = (visibleEnd - visibleStart + 1) / mapLength * 100;
        const isSelected = state.selectedFeatureIndex === index;
        const isVisible = visibleFeatureIndexes.has(index);
        if (app.isCdsFeature(feature)) {
          const lane = app.assignOverlapLane(cdsLaneEnds, visibleStart, visibleEnd);
          const marker = buildMarker(feature, index, "rgba(14, 165, 233, 0.86)", "pill");
          marker.className = `cursor-pointer ${isSelected ? "ring-2 ring-primary ring-offset-1 ring-offset-base-200 rounded" : ""} ${isVisible ? "opacity-100" : "opacity-25"}`;
          marker.style.left = `${left}%`;
          marker.style.width = `${Math.max(0.1, Math.min(width, 100 - left))}%`;
          marker.style.minWidth = "10px";
          marker.style.top = `${lane * (laneHeight + laneGap) + 2}px`;
          marker.style.height = `${laneHeight}px`;
          els.cdsFeatureTrack.appendChild(marker);
          cdsCount += 1;
          return;
        }

        if (!app.isPrimerOrMiscFeature(feature)) {
          return;
        }

        const isPrimer = app.isPrimerBindingFeature(feature);
        const isForward = app.normalizedStrandValue(feature.strand) !== -1;
        const lane = app.assignOverlapLane(primerMiscLaneEnds, visibleStart, visibleEnd);
        const color = isPrimer
          ? (isForward ? "rgba(37, 99, 235, 0.92)" : "rgba(220, 38, 38, 0.92)")
          : "rgba(16, 185, 129, 0.86)";
        const shape = isPrimer ? (isForward ? "arrow-forward" : "arrow-reverse") : "pill";
        const marker = buildMarker(feature, index, color, shape);
        marker.className = `cursor-pointer ${isSelected ? "ring-2 ring-primary ring-offset-1 ring-offset-base-200 rounded" : ""} ${isVisible ? "opacity-100" : "opacity-30"}`;
        marker.style.left = `${left}%`;
        marker.style.width = `${Math.max(0.1, Math.min(width, 100 - left))}%`;
        marker.style.minWidth = isPrimer ? "12px" : "10px";
        marker.style.top = `${lane * (laneHeight + laneGap) + 2}px`;
        marker.style.height = `${laneHeight}px`;
        els.primerMiscFeatureTrack.appendChild(marker);
        primerMiscCount += 1;
      });

      const cdsLaneCount = Math.max(1, cdsLaneEnds.length);
      const primerMiscLaneCount = Math.max(1, primerMiscLaneEnds.length);
      els.cdsFeatureTrackContainer.style.height = `${cdsLaneCount * (laneHeight + laneGap) + 6}px`;
      els.featureTrackContainer.style.height = `${primerMiscLaneCount * (laneHeight + laneGap) + 6}px`;
      els.primerCountSummary.textContent = `Map lanes: CDS ${cdsCount.toLocaleString()} across ${cdsLaneCount} lane(s) | primer_bind/misc_feature ${primerMiscCount.toLocaleString()} across ${primerMiscLaneCount} lane(s)`;
    };

    app.renderWindowFeatureOverlay = function renderWindowFeatureOverlay(record, start, end) {
      els.windowFeatureOverlay.innerHTML = "";
      els.windowFeatureOverlay.classList.add("window-feature-overlay-grid");

      if (!state.showWindowFeatureOverlay) {
        els.windowFeatureOverlaySection.classList.add("hidden");
        els.toggleWindowFeatureOverlayBtn.textContent = "Show feature lines";
        app.setWindowFeatureHoverText(null);
        return;
      }

      els.windowFeatureOverlaySection.classList.remove("hidden");
      els.toggleWindowFeatureOverlayBtn.textContent = "Hide feature lines";

      const windowLength = Math.max(1, end - start + 1);
      const visibleFeatureIndexes = new Set(app.getFilteredFeatureIndexes(record));
      const overlappingFeatures = record.features
        .map((feature, index) => ({ feature, index }))
        .filter(({ feature }) => {
          const bounds = app.featureBounds(feature);
          return (app.isCdsFeature(feature) || app.isPrimerOrMiscFeature(feature)) && !(bounds.end < start || bounds.start > end);
        })
        .sort((a, b) => {
          const aBounds = app.featureBounds(a.feature);
          const bBounds = app.featureBounds(b.feature);
          const startDelta = aBounds.start - bBounds.start;
          if (startDelta !== 0) {
            return startDelta;
          }
          return (aBounds.end - aBounds.start) - (bBounds.end - bBounds.start);
        });

      els.windowFeatureOverlayCount.textContent = `${overlappingFeatures.length.toLocaleString()} in view`;
      if (!overlappingFeatures.length) {
        els.windowFeatureOverlay.style.height = "2rem";
        const emptyState = document.createElement("div");
        emptyState.className = "absolute inset-0 flex items-center justify-center text-[11px] opacity-60";
        emptyState.textContent = "No CDS, primer_bind, or misc_feature annotations overlap the current sequence window.";
        els.windowFeatureOverlay.appendChild(emptyState);
        app.setWindowFeatureHoverText(null);
        return;
      }

      const cdsLaneEnds = [];
      const primerMiscLaneEnds = [];
      const overlayLaneHeight = 10;
      const overlayLaneGap = 4;
      const primerMiscLaneOffset = 5;

      overlappingFeatures.forEach(({ feature, index }) => {
        const bounds = app.featureBounds(feature);
        const visibleStart = Math.max(bounds.start, start);
        const visibleEnd = Math.min(bounds.end, end);
        const isCds = app.isCdsFeature(feature);
        const laneEnds = isCds ? cdsLaneEnds : primerMiscLaneEnds;
        const lane = app.assignOverlapLane(laneEnds, visibleStart, visibleEnd);

        const left = ((visibleStart - start) / windowLength) * 100;
        const width = ((visibleEnd - visibleStart + 1) / windowLength) * 100;
        const selected = state.selectedFeatureIndex === index;
        const visibleInFilter = visibleFeatureIndexes.has(index);
        const color = isCds
          ? "rgba(14, 165, 233, 0.92)"
          : (app.isPrimerBindingFeature(feature)
            ? (app.normalizedStrandValue(feature.strand) === -1 ? "rgba(220, 38, 38, 0.95)" : "rgba(37, 99, 235, 0.95)")
            : "rgba(16, 185, 129, 0.92)");

        const marker = document.createElement("button");
        marker.type = "button";
        marker.className = "absolute rounded-full";
        marker.style.left = `${left}%`;
        marker.style.width = `${Math.max(0.6, Math.min(width, 100 - left))}%`;
        marker.style.minWidth = "6px";
        marker.style.top = `${(isCds ? 8 : primerMiscLaneOffset + 8 + (Math.max(1, cdsLaneEnds.length) * (overlayLaneHeight + overlayLaneGap))) + lane * (overlayLaneHeight + overlayLaneGap)}px`;
        marker.style.height = `${overlayLaneHeight}px`;
        marker.style.backgroundColor = color;
        marker.style.border = "none";
        marker.style.padding = "0";
        marker.style.opacity = visibleInFilter ? "1" : "0.35";
        marker.style.boxShadow = selected
          ? "0 0 0 2px rgba(99, 102, 241, 0.35)"
          : "0 0 0 1px rgba(15, 23, 42, 0.18)";
        marker.style.cursor = "pointer";
        marker.title = `${feature.label} | ${app.featureDescription(feature)}`;
        marker.addEventListener("mouseenter", () => app.setWindowFeatureHoverText(feature));
        marker.addEventListener("focus", () => app.setWindowFeatureHoverText(feature));
        marker.addEventListener("mouseleave", () => {
          const selectedFeature = state.selectedFeatureIndex !== null
            ? record.features[state.selectedFeatureIndex]
            : null;
          app.setWindowFeatureHoverText(selectedFeature || null);
        });
        marker.addEventListener("blur", () => {
          const selectedFeature = state.selectedFeatureIndex !== null
            ? record.features[state.selectedFeatureIndex]
            : null;
          app.setWindowFeatureHoverText(selectedFeature || null);
        });
        marker.addEventListener("click", (event) => app.applyFeatureSelection(record, index, feature, event.shiftKey));
        els.windowFeatureOverlay.appendChild(marker);
      });

      const overlayHeight =
        16 +
        (Math.max(1, cdsLaneEnds.length) * (overlayLaneHeight + overlayLaneGap)) +
        primerMiscLaneOffset +
        (Math.max(1, primerMiscLaneEnds.length) * (overlayLaneHeight + overlayLaneGap));
      els.windowFeatureOverlay.style.height = `${overlayHeight}px`;
      const selectedFeature = state.selectedFeatureIndex !== null
        ? record.features[state.selectedFeatureIndex]
        : null;
      app.setWindowFeatureHoverText(selectedFeature || null);
    };

    app.renderMapAxis = function renderMapAxis(record) {
      els.mapTickTrack.innerHTML = "";
      const mapStart = state.mapStart;
      const mapEnd = Math.min(record.length, mapStart + state.mapWindowSize - 1);
      const mapLength = mapEnd - mapStart + 1;
      const step = Math.max(1, Math.floor(mapLength / 6));

      const tickValues = [];
      for (let value = mapStart; value <= mapEnd; value += step) {
        tickValues.push(value);
      }
      if (tickValues[tickValues.length - 1] !== mapEnd) {
        tickValues.push(mapEnd);
      }

      tickValues.forEach((value) => {
        const pct = ((value - mapStart) / mapLength) * 100;
        const tick = document.createElement("div");
        tick.className = "absolute top-0 h-3 w-px bg-base-content/35";
        tick.style.left = `${Math.min(100, Math.max(0, pct))}%`;
        els.mapTickTrack.appendChild(tick);

        const label = document.createElement("div");
        label.className = "absolute top-0 text-[10px] opacity-60";
        label.style.left = `${Math.min(100, Math.max(0, pct))}%`;
        label.style.transform = "translateX(-50%)";
        label.textContent = value.toLocaleString();
        els.mapTickTrack.appendChild(label);
      });

      els.mapLeftLabel.textContent = `${mapStart.toLocaleString()} bp`;
      els.mapRightLabel.textContent = `${mapEnd.toLocaleString()} bp`;
    };

    app.renderMapSelectionSummary = function renderMapSelectionSummary(record) {
      if (state.selectedFeatureIndex === null || !record.features[state.selectedFeatureIndex]) {
        els.mapSelectionSummary.textContent = "No feature selected on map.";
        els.mapSelectionActions.classList.add("hidden");
        return;
      }
      const feature = record.features[state.selectedFeatureIndex];
      els.mapSelectionSummary.textContent = `Selected: ${feature.label} (${feature.type}) ${feature.start.toLocaleString()}-${feature.end.toLocaleString()} | ${app.featureLength(feature).toLocaleString()} bp | ${app.strandLabel(feature.strand)}`;
      const canRemovePrimer = app.isPrimerBindingFeature(feature) && feature.source === "user" && Number(feature.feature_id) > 0;
      if (!canRemovePrimer) {
        els.mapSelectionActions.classList.add("hidden");
        return;
      }
      els.mapSelectionActions.classList.remove("hidden");
      els.removePrimerFeatureBtn.disabled = state.primerDeleteSubmitting;
      els.deletePrimerEverywhereBtn.disabled = state.primerDeleteSubmitting || !Number(feature.primer_id);
    };

    app.renderFeatureViewer = function renderFeatureViewer(record) {
      const tableIndexes = app.getTableFeatureIndexes(record);
      const selectedCount = state.mapSelectedFeatureIndexes.size;
      const modeLabel = state.showMapSelectedOnly ? "map-selected mode" : "all features mode";
      els.featureCount.textContent = `${tableIndexes.length.toLocaleString()} of ${record.features.length.toLocaleString()} features shown | map-selected: ${selectedCount.toLocaleString()} | ${modeLabel}`;
      if (!tableIndexes.length) {
        const emptyMessage = state.showMapSelectedOnly
          ? "No map-selected features match the current filter."
          : "No features match the current filter.";
        els.featureTableBody.innerHTML = `<tr><td colspan="6" class="text-center opacity-70 py-3">${emptyMessage}</td></tr>`;
        els.featurePageInfo.textContent = "Page 0 of 0";
        els.featurePrevPageBtn.disabled = true;
        els.featureNextPageBtn.disabled = true;
        return;
      }

      if (state.selectedFeatureIndex === null || !tableIndexes.includes(state.selectedFeatureIndex)) {
        state.selectedFeatureIndex = tableIndexes[0];
      }

      if (state.pendingFocusFeatureIndex !== null && tableIndexes.includes(state.pendingFocusFeatureIndex)) {
        state.selectedFeatureIndex = state.pendingFocusFeatureIndex;
        const focusedPosition = tableIndexes.indexOf(state.pendingFocusFeatureIndex);
        state.tablePage = Math.floor(focusedPosition / state.tablePageSize) + 1;
        state.pendingFocusFeatureIndex = null;
      } else if (state.pendingFocusFeatureIndex !== null && !tableIndexes.includes(state.pendingFocusFeatureIndex)) {
        state.pendingFocusFeatureIndex = null;
      }

      const totalPages = Math.max(1, Math.ceil(tableIndexes.length / state.tablePageSize));
      state.tablePage = Math.min(Math.max(1, state.tablePage), totalPages);
      const pageIndexes = tableIndexes.slice(
        (state.tablePage - 1) * state.tablePageSize,
        state.tablePage * state.tablePageSize
      );

      els.featurePageInfo.textContent = `Page ${state.tablePage.toLocaleString()} of ${totalPages.toLocaleString()}`;
      els.featurePrevPageBtn.disabled = state.tablePage <= 1;
      els.featureNextPageBtn.disabled = state.tablePage >= totalPages;

      els.featureTableBody.innerHTML = "";
      pageIndexes.forEach((index) => {
        const feature = record.features[index];
        const row = document.createElement("tr");
        row.className = `cursor-pointer ${state.selectedFeatureIndex === index ? "bg-base-200" : ""}`;
        row.innerHTML = `
          <td class="max-w-56 truncate" title="${feature.label}">${feature.label}</td>
          <td>${feature.type}</td>
          <td class="font-mono">${feature.start.toLocaleString()}</td>
          <td class="font-mono">${feature.end.toLocaleString()}</td>
          <td class="font-mono">${app.featureLength(feature).toLocaleString()} bp</td>
          <td>${app.strandLabel(feature.strand)}</td>
        `;
        row.addEventListener("click", () => {
          state.selectedFeatureIndex = index;
          state.pendingFocusFeatureIndex = index;
          app.jumpToFeature(record, feature);
          app.render();
        });
        els.featureTableBody.appendChild(row);
      });
    };

    app.displayMapSelectedFeaturesInTable = function displayMapSelectedFeaturesInTable() {
      state.showMapSelectedOnly = true;
      state.tablePage = 1;
      app.render();
    };
  };
})();
