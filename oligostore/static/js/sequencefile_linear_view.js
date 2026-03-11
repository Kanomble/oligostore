document.addEventListener("DOMContentLoaded", () => {
  const raw = document.getElementById("sequence-records-data");
  const config = document.getElementById("sequence-linear-config");
  if (!raw || !config) {
    return;
  }

    const recordSummaries = JSON.parse(raw.textContent);
    if (!recordSummaries.length) {
      return;
    }
    const recordDetails = new Map();
    const recordDataUrl = config.dataset.recordDataUrl;
    const createPrimerUrl = config.dataset.createPrimerUrl;
    const deletePrimerUrl = config.dataset.deletePrimerUrl;
    const analyzePrimerUrl = config.dataset.analyzePrimerUrl;

    const recordSelect = document.getElementById("recordSelect");
    const windowSizeSlider = document.getElementById("windowSizeSlider");
    const startSlider = document.getElementById("startSlider");
    const zoomInBtn = document.getElementById("zoomInBtn");
    const zoomOutBtn = document.getElementById("zoomOutBtn");
    const resetBtn = document.getElementById("resetBtn");
    const mapZoomInBtn = document.getElementById("mapZoomInBtn");
    const mapZoomOutBtn = document.getElementById("mapZoomOutBtn");
    const mapZoomResetBtn = document.getElementById("mapZoomResetBtn");
    const mapStartSlider = document.getElementById("mapStartSlider");
    const mapRangeLabel = document.getElementById("mapRangeLabel");
    const mapTickTrack = document.getElementById("mapTickTrack");
    const mapLeftLabel = document.getElementById("mapLeftLabel");
    const mapRightLabel = document.getElementById("mapRightLabel");
    const primerCountSummary = document.getElementById("primerCountSummary");
    const mapSelectionSummary = document.getElementById("mapSelectionSummary");
    const mapSelectionActions = document.getElementById("mapSelectionActions");
    const removePrimerFeatureBtn = document.getElementById("removePrimerFeatureBtn");
    const deletePrimerEverywhereBtn = document.getElementById("deletePrimerEverywhereBtn");
    const pcrProductSummary = document.getElementById("pcrProductSummary");
    const pcrProductSequence = document.getElementById("pcrProductSequence");
    const windowLabel = document.getElementById("windowLabel");
    const windowSizeHint = document.getElementById("windowSizeHint");
    const viewportBox = document.getElementById("viewportBox");
    const featureTrackContainer = document.getElementById("featureTrackContainer");
    const featureTrack = document.getElementById("featureTrack");
    const forwardPrimerTrack = document.getElementById("forwardPrimerTrack");
    const reversePrimerTrack = document.getElementById("reversePrimerTrack");
    const featureLegend = document.getElementById("featureLegend");
    const featureSearchInput = document.getElementById("featureSearchInput");
    const featureCount = document.getElementById("featureCount");
    const featureTableBody = document.getElementById("featureTableBody");
    const showMapSelectedBtn = document.getElementById("showMapSelectedBtn");
    const showAllFeaturesBtn = document.getElementById("showAllFeaturesBtn");
    const clearMapSelectionBtn = document.getElementById("clearMapSelectionBtn");
    const featurePageSizeSelect = document.getElementById("featurePageSizeSelect");
    const featurePageInfo = document.getElementById("featurePageInfo");
    const featurePrevPageBtn = document.getElementById("featurePrevPageBtn");
    const featureNextPageBtn = document.getElementById("featureNextPageBtn");
    const selectTopRestrictionEnzymesBtn = document.getElementById("selectTopRestrictionEnzymesBtn");
    const clearRestrictionSelectionBtn = document.getElementById("clearRestrictionSelectionBtn");
    const restrictionEnzymeSearchInput = document.getElementById("restrictionEnzymeSearchInput");
    const restrictionEnzymeTableBody = document.getElementById("restrictionEnzymeTableBody");
    const restrictionEnzymePageInfo = document.getElementById("restrictionEnzymePageInfo");
    const restrictionEnzymePrevPageBtn = document.getElementById("restrictionEnzymePrevPageBtn");
    const restrictionEnzymeNextPageBtn = document.getElementById("restrictionEnzymeNextPageBtn");
    const restrictionSelectionCount = document.getElementById("restrictionSelectionCount");
    const restrictionSiteSummary = document.getElementById("restrictionSiteSummary");
    const sequenceWindow = document.getElementById("sequenceWindow");
    const primerSelectionMenu = document.getElementById("primerSelectionMenu");
    const closePrimerSelectionMenuBtn = document.getElementById("closePrimerSelectionMenuBtn");
    const primerSelectionSummary = document.getElementById("primerSelectionSummary");
    const primerSelectionSequence = document.getElementById("primerSelectionSequence");
    const primerSelectionTm = document.getElementById("primerSelectionTm");
    const primerSelectionGc = document.getElementById("primerSelectionGc");
    const primerSelectionHairpin = document.getElementById("primerSelectionHairpin");
    const primerSelectionSelfDimer = document.getElementById("primerSelectionSelfDimer");
    const primerNameInput = document.getElementById("primerNameInput");
    const primerOverhangInput = document.getElementById("primerOverhangInput");
    const savePrimerToOligostoreCheckbox = document.getElementById("savePrimerToOligostoreCheckbox");
    const attachPrimerAsFeatureCheckbox = document.getElementById("attachPrimerAsFeatureCheckbox");
    const reverseComplementPrimerBtn = document.getElementById("reverseComplementPrimerBtn");
    const savePrimerFromWindowBtn = document.getElementById("savePrimerFromWindowBtn");
    const primerCreateStatus = document.getElementById("primerCreateStatus");
    const toggleWindowFeatureOverlayBtn = document.getElementById("toggleWindowFeatureOverlayBtn");
    const windowFeatureOverlaySection = document.getElementById("windowFeatureOverlaySection");
    const windowFeatureOverlay = document.getElementById("windowFeatureOverlay");
    const windowFeatureOverlayCount = document.getElementById("windowFeatureOverlayCount");
    const windowFeatureHoverDescription = document.getElementById("windowFeatureHoverDescription");

    const state = {
      recordIndex: 0,
      windowSize: 200,
      start: 1,
      mapWindowSize: 1000,
      mapStart: 1,
      selectedFeatureIndex: null,
      featureQuery: "",
      tablePage: 1,
      tablePageSize: 10,
      showMapSelectedOnly: false,
      mapSelectedFeatureIndexes: new Set(),
      pendingFocusFeatureIndex: null,
      selectedForwardPrimerIndex: null,
      selectedReversePrimerIndex: null,
      selectedRestrictionEnzymes: new Set(),
      restrictionEnzymeQuery: "",
      restrictionTablePage: 1,
      restrictionTablePageSize: 10,
      primerSubmitting: false,
      loadingRecord: false,
      loadError: "",
      pendingRegionRequest: null,
      selectedPrimerCandidate: null,
      primerAnalysisLoading: false,
      primerDeleteSubmitting: false,
      showWindowFeatureOverlay: true,
    };

    const MIN_WINDOW_BP = 1;

    function getRecord(index = state.recordIndex) {
      return recordDetails.get(index) || null;
    }

    function getRecordSummary() {
      return recordSummaries[state.recordIndex] || null;
    }

    function getCurrentRecordLength() {
      const loaded = getRecord();
      if (loaded) {
        return loaded.length;
      }
      const summary = getRecordSummary();
      return summary ? summary.length : 1;
    }

    function getActiveRegionBounds() {
      const recordLength = getCurrentRecordLength();
      const sequenceEnd = Math.min(recordLength, state.start + state.windowSize - 1);
      const mapEnd = Math.min(recordLength, state.mapStart + state.mapWindowSize - 1);
      const regionStart = Math.max(1, Math.min(state.start, state.mapStart));
      const regionEnd = Math.max(regionStart, Math.max(sequenceEnd, mapEnd));
      return { regionStart, regionEnd };
    }

    function isRangeCovered(record, neededStart, neededEnd) {
      if (!record) {
        return false;
      }
      return record.region_start <= neededStart && record.region_end >= neededEnd;
    }

    async function ensureRecordRegionLoaded(index, neededStart, neededEnd) {
      const cached = getRecord(index);
      if (isRangeCovered(cached, neededStart, neededEnd)) {
        return;
      }
      if (state.loadingRecord) {
        state.pendingRegionRequest = { index, start: neededStart, end: neededEnd };
        return;
      }
      if (state.recordIndex === index) {
        state.loadingRecord = true;
        state.loadError = "";
        render();
      }
      try {
        const response = await fetch(
          `${recordDataUrl}?record_index=${index}&start=${neededStart}&end=${neededEnd}`
        );
        if (!response.ok) {
          throw new Error(`HTTP ${response.status}`);
        }
        const payload = await response.json();
        recordDetails.set(index, payload);
      } catch (error) {
        if (state.recordIndex === index) {
          state.loadError = `Could not load sequence record data (${error.message}).`;
        }
      } finally {
        if (state.recordIndex === index) {
          state.loadingRecord = false;
          render();
        }
        const pending = state.pendingRegionRequest;
        state.pendingRegionRequest = null;
        if (pending) {
          void ensureRecordRegionLoaded(pending.index, pending.start, pending.end);
        }
      }
    }

    function clampWindowToRecord() {
      const recordLength = getCurrentRecordLength();
      const maxWindow = Math.min(5000, Math.max(MIN_WINDOW_BP, recordLength));
      state.windowSize = Math.min(Math.max(MIN_WINDOW_BP, state.windowSize), maxWindow);
      if (state.start > recordLength) {
        state.start = 1;
      }
      const maxStart = Math.max(1, recordLength - state.windowSize + 1);
      state.start = Math.min(state.start, maxStart);
      state.start = Math.max(1, state.start);
    }

    function defaultMapWindowSize(recordLength) {
      return Math.max(50, Math.min(5000, Math.max(1, recordLength)));
    }

    function clampMapToRecord() {
      const recordLength = getCurrentRecordLength();
      const maxMapWindow = Math.max(50, recordLength);
      state.mapWindowSize = Math.min(Math.max(50, state.mapWindowSize), maxMapWindow);
      if (state.mapStart > recordLength) {
        state.mapStart = 1;
      }
      const maxMapStart = Math.max(1, recordLength - state.mapWindowSize + 1);
      state.mapStart = Math.min(state.mapStart, maxMapStart);
      state.mapStart = Math.max(1, state.mapStart);
    }

    function formatFeatures(features) {
      if (!features.length) {
        return "No annotated features in this record.";
      }
      const byType = {};
      for (const feature of features) {
        byType[feature.type] = (byType[feature.type] || 0) + 1;
      }
      return Object.entries(byType)
        .sort((a, b) => a[0].localeCompare(b[0]))
        .map(([type, count]) => `${type}: ${count}`)
        .join(" | ");
    }

    function baseSpan(base, highlightRestriction = false, position = null, strand = "forward") {
      const value = base.toUpperCase();
      const restrictionClass = highlightRestriction ? " restriction-hit" : "";
      const metadata = Number.isFinite(position)
        ? ` data-sequence-base="1" data-position="${position}" data-strand="${strand}"`
        : "";
      const classPrefix = Number.isFinite(position) ? "sequence-base " : "";
      if (value === "A") return `<span${metadata} class="${classPrefix}base-a${restrictionClass}">A</span>`;
      if (value === "C") return `<span${metadata} class="${classPrefix}base-c${restrictionClass}">C</span>`;
      if (value === "G") return `<span${metadata} class="${classPrefix}base-g${restrictionClass}">G</span>`;
      if (value === "T") return `<span${metadata} class="${classPrefix}base-t${restrictionClass}">T</span>`;
      return `<span${metadata} class="${classPrefix}base-n${restrictionClass}">${value}</span>`;
    }

    const COMPLEMENT_BY_BASE = {
      A: "T",
      T: "A",
      C: "G",
      G: "C",
    };
    const SEQUENCE_WINDOW_CHUNK_SIZE = 10;
    const SEQUENCE_WINDOW_LINE_SIZE = 100;

    function complementBase(base) {
      const value = String(base || "").toUpperCase();
      return COMPLEMENT_BY_BASE[value] || "N";
    }

    function escapeHtml(value) {
      return String(value)
        .replace(/&/g, "&amp;")
        .replace(/</g, "&lt;")
        .replace(/>/g, "&gt;")
        .replace(/"/g, "&quot;")
        .replace(/'/g, "&#39;");
    }

    function getRestrictionHitsInWindow(record, start, end) {
      return (record.restriction_sites || []).filter((site) => !(site.end < start || site.start > end));
    }

    function getRestrictionEnzymeStats(record) {
      const stats = {};
      (record.restriction_sites || []).forEach((site) => {
        if (!stats[site.enzyme]) {
          stats[site.enzyme] = {
            enzyme: site.enzyme,
            count: 0,
            site: site.site || "",
            cutOffset: Number(site.cut_offset),
          };
        }
        stats[site.enzyme].count += 1;
        if (!stats[site.enzyme].site && site.site) {
          stats[site.enzyme].site = site.site;
        }
        if (!Number.isFinite(stats[site.enzyme].cutOffset) && Number.isFinite(Number(site.cut_offset))) {
          stats[site.enzyme].cutOffset = Number(site.cut_offset);
        }
      });
      return stats;
    }

    function formatNRun(count) {
      if (count <= 0) {
        return "";
      }
      if (count <= 8) {
        return "N".repeat(count);
      }
      return `N(${count})`;
    }

    function formatCutSite(site, cutOffset) {
      const motif = String(site || "").toUpperCase();
      if (!motif) {
        return "-";
      }
      if (!Number.isFinite(cutOffset)) {
        return motif;
      }
      const offset = Math.trunc(cutOffset);
      if (offset >= 0 && offset <= motif.length) {
        return `${motif.slice(0, offset)}/${motif.slice(offset)}`;
      }
      if (offset < 0) {
        return `${formatNRun(Math.abs(offset))}/${motif}`;
      }
      return `${motif}/${formatNRun(offset - motif.length)}`;
    }

    function getSortedRestrictionEnzymes(record) {
      const stats = getRestrictionEnzymeStats(record);
      return Object.values(stats)
        .sort((a, b) => {
          const countDelta = b.count - a.count;
          if (countDelta !== 0) {
            return countDelta;
          }
          return a.enzyme.localeCompare(b.enzyme);
        });
    }

    function applyTopRestrictionSelection(record, limit = 5) {
      const sorted = getSortedRestrictionEnzymes(record);
      state.selectedRestrictionEnzymes = new Set(
        sorted.slice(0, limit).map((item) => item.enzyme)
      );
    }

    function renderRestrictionEnzymeTable(record) {
      const sorted = getSortedRestrictionEnzymes(record);
      const available = new Set(sorted.map((item) => item.enzyme));
      state.selectedRestrictionEnzymes = new Set(
        [...state.selectedRestrictionEnzymes].filter((enzyme) => available.has(enzyme))
      );
      const query = normalize(state.restrictionEnzymeQuery).trim();
      const selectedItems = sorted.filter((item) => state.selectedRestrictionEnzymes.has(item.enzyme));
      const unselectedItems = sorted.filter((item) => !state.selectedRestrictionEnzymes.has(item.enzyme));
      const filteredUnselectedItems = !query
        ? unselectedItems
        : unselectedItems.filter(({ enzyme }) => normalize(enzyme).includes(query));

      if (!sorted.length) {
        restrictionEnzymeTableBody.innerHTML = '<tr><td colspan="4" class="text-center opacity-70 py-3">No restriction sites in this record.</td></tr>';
        restrictionSelectionCount.textContent = "0 selected";
        restrictionEnzymePageInfo.textContent = "Page 0 of 0";
        restrictionEnzymePrevPageBtn.disabled = true;
        restrictionEnzymeNextPageBtn.disabled = true;
        return;
      }

      if (!selectedItems.length && !filteredUnselectedItems.length) {
        restrictionEnzymeTableBody.innerHTML = '<tr><td colspan="4" class="text-center opacity-70 py-3">No enzymes match this search.</td></tr>';
        restrictionSelectionCount.textContent = `${state.selectedRestrictionEnzymes.size.toLocaleString()} selected | 0 matching`;
        restrictionEnzymePageInfo.textContent = "Page 0 of 0";
        restrictionEnzymePrevPageBtn.disabled = true;
        restrictionEnzymeNextPageBtn.disabled = true;
        return;
      }

      const totalPages = Math.max(1, Math.ceil(filteredUnselectedItems.length / state.restrictionTablePageSize));
      state.restrictionTablePage = Math.min(Math.max(1, state.restrictionTablePage), totalPages);
      const pageStart = (state.restrictionTablePage - 1) * state.restrictionTablePageSize;
      const pageEnd = pageStart + state.restrictionTablePageSize;
      const pageItems = filteredUnselectedItems.slice(pageStart, pageEnd);
      const rowsToRender = [...selectedItems, ...pageItems];

      restrictionEnzymeTableBody.innerHTML = "";
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
          render();
        });

        const selectCell = document.createElement("td");
        selectCell.appendChild(input);
        const enzymeCell = document.createElement("td");
        enzymeCell.className = "font-mono";
        enzymeCell.textContent = enzyme;
        const siteCell = document.createElement("td");
        siteCell.className = "font-mono";
        const cutSiteLabel = formatCutSite(site, cutOffset);
        siteCell.textContent = cutSiteLabel;
        siteCell.title = `${site || "-"} | cut offset: ${Number.isFinite(cutOffset) ? cutOffset : "n/a"}`;
        const countCell = document.createElement("td");
        countCell.textContent = count.toLocaleString();

        row.appendChild(selectCell);
        row.appendChild(enzymeCell);
        row.appendChild(siteCell);
        row.appendChild(countCell);
        restrictionEnzymeTableBody.appendChild(row);
      });

      const matchingTotal = selectedItems.length + filteredUnselectedItems.length;
      restrictionSelectionCount.textContent = `${state.selectedRestrictionEnzymes.size.toLocaleString()} selected | ${matchingTotal.toLocaleString()} matching`;
      restrictionEnzymePageInfo.textContent = `Page ${state.restrictionTablePage.toLocaleString()} of ${totalPages.toLocaleString()}`;
      restrictionEnzymePrevPageBtn.disabled = state.restrictionTablePage <= 1;
      restrictionEnzymeNextPageBtn.disabled = state.restrictionTablePage >= totalPages;
    }

    function normalize(value) {
      return String(value || "").toLowerCase();
    }

    function getFilteredFeatureIndexes(record) {
      const query = normalize(state.featureQuery).trim();
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
          ].map(normalize).join(" ");
          return haystack.includes(query);
        })
        .map(({ index }) => index);
    }

    function getTableFeatureIndexes(record) {
      const filtered = getFilteredFeatureIndexes(record);
      if (!state.showMapSelectedOnly) {
        return filtered;
      }
      return filtered.filter((index) => state.mapSelectedFeatureIndexes.has(index));
    }

    function featureLength(feature) {
      const bounds = featureBounds(feature);
      return bounds.end - bounds.start + 1;
    }

    function featureBounds(feature) {
      const rawStart = Number(feature.start);
      const rawEnd = Number(feature.end);
      const safeStart = Number.isFinite(rawStart) ? rawStart : 1;
      const safeEnd = Number.isFinite(rawEnd) ? rawEnd : safeStart;
      return {
        start: Math.min(safeStart, safeEnd),
        end: Math.max(safeStart, safeEnd),
      };
    }

    function featureDescription(feature) {
      const description = String(
        feature.description || feature.note || feature.product || feature.label || feature.type || ""
      ).trim();
      return description || "No description available.";
    }

    function isPrimerBindingFeature(feature) {
      const typeValue = normalize(feature.type);
      const labelValue = normalize(feature.label);
      const noteValue = normalize(feature.note);
      if (typeValue.includes("primer") || typeValue.includes("primer_bind")) {
        return true;
      }
      if (labelValue.includes("primer") || noteValue.includes("primer")) {
        return true;
      }
      if (labelValue.includes("_fw") || labelValue.includes("_rv") || labelValue.includes("_rev")) {
        return true;
      }
      return false;
    }

    function normalizedStrandValue(strand) {
      if (strand === -1 || strand === 1) {
        return strand;
      }
      const numeric = Number(strand);
      if (numeric === -1 || numeric === 1) {
        return numeric;
      }
      const strandText = normalize(strand);
      if (strandText.includes("reverse") || strandText.includes("complement") || strandText === "-") {
        return -1;
      }
      if (strandText.includes("forward") || strandText === "+") {
        return 1;
      }
      return 1;
    }

    function strandLabel(strand) {
      if (strand === 1) return "Forward (+)";
      if (strand === -1) return "Reverse (-)";
      return "Unknown";
    }

    function jumpToFeature(record, feature) {
      const length = featureLength(feature);
      const bounds = featureBounds(feature);
      const maxWindow = Math.min(5000, record.length);
      state.windowSize = Math.max(MIN_WINDOW_BP, Math.min(maxWindow, length));
      const midpoint = Math.floor((bounds.start + bounds.end) / 2);
      const centeredStart = midpoint - Math.floor(state.windowSize / 2);
      const maxStart = Math.max(1, record.length - state.windowSize + 1);
      state.start = Math.max(1, Math.min(centeredStart, maxStart));
      const mapMidpoint = Math.floor((bounds.start + bounds.end) / 2);
      const centeredMapStart = mapMidpoint - Math.floor(state.mapWindowSize / 2);
      const maxMapStart = Math.max(1, record.length - state.mapWindowSize + 1);
      state.mapStart = Math.max(1, Math.min(centeredMapStart, maxMapStart));
    }

    function applyFeatureSelection(record, index, feature, isShift) {
      state.selectedFeatureIndex = index;
      if (isShift) {
        state.mapSelectedFeatureIndexes.add(index);
      } else {
        state.mapSelectedFeatureIndexes = new Set([index]);
      }
      state.showMapSelectedOnly = true;
      state.pendingFocusFeatureIndex = index;
      if (isShift && isPrimerBindingFeature(feature)) {
        if (state.selectedForwardPrimerIndex !== null && state.selectedReversePrimerIndex !== null) {
          state.selectedForwardPrimerIndex = null;
          state.selectedReversePrimerIndex = null;
          state.mapSelectedFeatureIndexes = new Set([index]);
        }
        const strand = normalizedStrandValue(feature.strand);
        if (strand === -1) {
          state.selectedReversePrimerIndex = index;
        } else {
          state.selectedForwardPrimerIndex = index;
        }
      }
      state.tablePage = 1;
      const didFocusPCRProduct = isShift && isPrimerBindingFeature(feature) && focusWindowOnPCRProduct(record);
      if (!didFocusPCRProduct) {
        jumpToFeature(record, feature);
      }
      render();
    }

    function setWindowFeatureHoverText(feature) {
      if (!feature) {
        windowFeatureHoverDescription.textContent = "Hover a feature line to inspect its description.";
        return;
      }
      windowFeatureHoverDescription.textContent =
        `${feature.label} (${feature.type}) ${feature.start.toLocaleString()}-${feature.end.toLocaleString()} | ` +
        `${featureLength(feature).toLocaleString()} bp | ${strandLabel(feature.strand)} | ${featureDescription(feature)}`;
    }

    function computePCRProduct(record) {
      if (state.selectedForwardPrimerIndex === null || state.selectedReversePrimerIndex === null) {
        return null;
      }
      const forwardPrimer = record.features[state.selectedForwardPrimerIndex];
      const reversePrimer = record.features[state.selectedReversePrimerIndex];
      if (!forwardPrimer || !reversePrimer) {
        return null;
      }
      const fwdBounds = featureBounds(forwardPrimer);
      const revBounds = featureBounds(reversePrimer);
      const productStart = fwdBounds.start;
      const productEnd = revBounds.end;
      if (productStart > productEnd) {
        return {
          valid: false,
          reason: "Forward primer is downstream of reverse primer in this linear view.",
        };
      }
      const sequence = getSequenceSlice(record, productStart, productEnd);
      if (sequence === null) {
        return {
          valid: false,
          reason: "PCR product sequence is outside the loaded region. Pan/zoom to load that region.",
        };
      }
      return {
        valid: true,
        start: productStart,
        end: productEnd,
        length: sequence.length,
        sequence,
        forwardLabel: forwardPrimer.label,
        reverseLabel: reversePrimer.label,
      };
    }

    function renderPCRProduct(record) {
      const product = computePCRProduct(record);
      if (!product) {
        pcrProductSummary.textContent = "PCR product: shift-click one forward and one reverse primer to generate a candidate amplicon.";
        pcrProductSequence.textContent = "";
        return;
      }
      if (!product.valid) {
        pcrProductSummary.textContent = `PCR product unavailable: ${product.reason}`;
        pcrProductSequence.textContent = "";
        return;
      }
      pcrProductSummary.textContent = `PCR product: ${product.start.toLocaleString()}-${product.end.toLocaleString()} (${product.length.toLocaleString()} bp) | Fwd: ${product.forwardLabel} | Rev: ${product.reverseLabel}`;
      pcrProductSequence.textContent = product.sequence;
    }

    function focusWindowOnPCRProduct(record) {
      const product = computePCRProduct(record);
      if (!product || !product.valid) {
        return false;
      }
      const maxWindow = Math.min(5000, record.length);
      state.windowSize = Math.max(MIN_WINDOW_BP, Math.min(maxWindow, product.length));
      state.start = Math.max(1, Math.min(product.start, record.length - state.windowSize + 1));
      const midpoint = Math.floor((product.start + product.end) / 2);
      const centeredMapStart = midpoint - Math.floor(state.mapWindowSize / 2);
      const maxMapStart = Math.max(1, record.length - state.mapWindowSize + 1);
      state.mapStart = Math.max(1, Math.min(centeredMapStart, maxMapStart));
      return true;
    }

    function displayMapSelectedFeaturesInTable() {
      state.showMapSelectedOnly = true;
      state.tablePage = 1;
      render();
    }

    function getHighlightedRestrictionPositions(windowRestrictionHits, start, end) {
      const highlightedPositions = new Set();
      windowRestrictionHits.forEach((site) => {
        const siteStart = Math.max(start, site.start);
        const siteEnd = Math.min(end, site.end);
        for (let position = siteStart; position <= siteEnd; position += 1) {
          highlightedPositions.add(position);
        }
      });
      return highlightedPositions;
    }

    function renderStrandChunks(lineSegment, lineStart, highlightedPositions, strand, baseTransform = (base) => base) {
      const chunks = [];
      for (let offset = 0; offset < lineSegment.length; offset += SEQUENCE_WINDOW_CHUNK_SIZE) {
        const chunk = lineSegment.slice(offset, offset + SEQUENCE_WINDOW_CHUNK_SIZE);
        const chunkBaseIndex = lineStart + offset;
        const renderedChunk = chunk.split("").map((base, chunkOffset) => {
          const position = chunkBaseIndex + chunkOffset;
          return baseSpan(baseTransform(base), highlightedPositions.has(position), position, strand);
        }).join("");
        chunks.push(renderedChunk);
      }
      return chunks.join(" ");
    }

    function renderLineRestrictionSites(windowRestrictionHits, lineStart, lineEnd) {
      return windowRestrictionHits
        .filter((site) => !(site.end < lineStart || site.start > lineEnd))
        .map((site) => `${escapeHtml(site.enzyme)} ${site.start.toLocaleString()}-${site.end.toLocaleString()}`);
    }

    function sequenceOffsetToDisplayOffset(offset) {
      return offset + Math.floor(offset / SEQUENCE_WINDOW_CHUNK_SIZE);
    }

    function sequenceDisplayLength(lineLength) {
      if (lineLength <= 0) {
        return 0;
      }
      return lineLength + Math.floor((lineLength - 1) / SEQUENCE_WINDOW_CHUNK_SIZE);
    }

    function buildRestrictionAnnotationRows(windowRestrictionHits, lineStart, lineEnd) {
      const lineLength = lineEnd - lineStart + 1;
      const displayLength = sequenceDisplayLength(lineLength);
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
        const markerIndex = sequenceOffsetToDisplayOffset(cutOffsetClamped);
        markerChars[markerIndex] = "|";

        const label = String(site.enzyme || "").trim();
        if (!label) {
          return;
        }
        const labelStart = sequenceOffsetToDisplayOffset(startOffset);
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
    }

    function getLineNumberLayout(end) {
      const width = String(Number(end).toLocaleString()).length;
      const prefix = " ".repeat(width + 2);
      return { width, prefix };
    }

    function getCsrfToken() {
      const match = document.cookie.match(/(?:^|; )csrftoken=([^;]+)/);
      return match ? decodeURIComponent(match[1]) : "";
    }

    function reverseComplementSequence(sequence) {
      return String(sequence || "")
        .toUpperCase()
        .split("")
        .reverse()
        .map((base) => complementBase(base))
        .join("");
    }

    function setPrimerStatus(message, isError = false) {
      primerCreateStatus.textContent = message;
      primerCreateStatus.className = `text-xs mt-2 ${isError ? "text-error" : "text-success"}`;
    }

    function getCurrentRecordId() {
      const summary = getRecordSummary();
      return summary ? String(summary.id) : "";
    }

    function closePrimerSelectionMenu() {
      primerSelectionMenu.classList.add("hidden");
      state.selectedPrimerCandidate = null;
      state.primerAnalysisLoading = false;
      savePrimerFromWindowBtn.disabled = true;
      setPrimerStatus("");
    }

    function setPrimerAnalysisPlaceholders(message = "-") {
      primerSelectionTm.textContent = message;
      primerSelectionGc.textContent = message;
      primerSelectionHairpin.textContent = message;
      primerSelectionSelfDimer.textContent = message;
    }

    function formatPrimerSelectionName(candidate) {
      const record = getRecord();
      const safeRecord = String((record && record.id) || `record_${state.recordIndex + 1}`).replace(/[^A-Za-z0-9_-]/g, "_");
      const direction = candidate.strand === -1 ? "R" : "F";
      return `${safeRecord}_${candidate.start}_${candidate.end}_${direction}`;
    }

    function getSelectedPrimerCandidate() {
      const selection = window.getSelection();
      if (!selection || selection.rangeCount === 0 || selection.isCollapsed) {
        return null;
      }
      const range = selection.getRangeAt(0);
      if (!sequenceWindow.contains(range.commonAncestorContainer)) {
        return null;
      }

      const selectedBaseSpans = Array.from(
        sequenceWindow.querySelectorAll('span[data-sequence-base="1"]')
      ).filter((span) => {
        try {
          return range.intersectsNode(span);
        } catch (error) {
          return false;
        }
      });

      if (!selectedBaseSpans.length) {
        return null;
      }

      const strandValues = [...new Set(selectedBaseSpans.map((span) => span.dataset.strand))];
      if (strandValues.length !== 1) {
        return { error: "Select bases from only one strand before right-clicking." };
      }

      const strandName = strandValues[0];
      const strand = strandName === "reverse" ? -1 : 1;
      const positions = selectedBaseSpans.map((span) => Number(span.dataset.position));
      if (positions.some((value) => !Number.isFinite(value))) {
        return { error: "Could not determine the selected base positions." };
      }

      const uniquePositions = [...new Set(positions)].sort((a, b) => a - b);
      if (!uniquePositions.length) {
        return null;
      }
      const expectedLength = uniquePositions[uniquePositions.length - 1] - uniquePositions[0] + 1;
      if (expectedLength !== uniquePositions.length) {
        return { error: "Select one continuous primer region." };
      }

      const displayedSequence = selectedBaseSpans
        .map((span) => String(span.textContent || "").trim().toUpperCase())
        .join("");
      if (!displayedSequence || displayedSequence.length !== uniquePositions.length) {
        return { error: "Select complete bases inside the sequence window." };
      }

      return {
        start: uniquePositions[0],
        end: uniquePositions[uniquePositions.length - 1],
        strand,
        displaySequence: displayedSequence,
        sequence: strand === -1
          ? displayedSequence.split("").reverse().join("")
          : displayedSequence,
      };
    }

    async function analyzePrimerCandidate(candidate) {
      state.primerAnalysisLoading = true;
      setPrimerAnalysisPlaceholders("Loading...");
      try {
        const response = await fetch(analyzePrimerUrl, {
          method: "POST",
          headers: {
            "Content-Type": "application/x-www-form-urlencoded;charset=UTF-8",
            "X-CSRFToken": getCsrfToken(),
          },
          body: new URLSearchParams({
            sequence: candidate.sequence,
          }),
        });
        const payload = await response.json().catch(() => ({}));
        if (!response.ok) {
          throw new Error(payload.error || `HTTP ${response.status}`);
        }
        primerSelectionTm.textContent = `${payload.tm} deg C`;
        primerSelectionGc.textContent = `${Math.round(Number(payload.gc_content || 0) * 100)}%`;
        primerSelectionHairpin.textContent = payload.hairpin
          ? `Possible (${payload.hairpin_dg} kcal/mol)`
          : `Not detected (${payload.hairpin_dg} kcal/mol)`;
        primerSelectionSelfDimer.textContent = payload.self_dimer
          ? `Possible (${payload.self_dimer_dg} kcal/mol)`
          : `Not detected (${payload.self_dimer_dg} kcal/mol)`;
      } catch (error) {
        setPrimerAnalysisPlaceholders("Unavailable");
        setPrimerStatus(`Could not analyze selection: ${error.message}`, true);
      } finally {
        state.primerAnalysisLoading = false;
      }
    }

    function openPrimerSelectionMenu(candidate, x, y) {
      state.selectedPrimerCandidate = candidate;
      primerSelectionSummary.textContent = `${candidate.start.toLocaleString()}-${candidate.end.toLocaleString()} | ${candidate.sequence.length.toLocaleString()} bp | ${strandLabel(candidate.strand)}`;
      primerSelectionSequence.textContent = candidate.sequence;
      primerNameInput.value = formatPrimerSelectionName(candidate);
      primerOverhangInput.value = "";
      savePrimerToOligostoreCheckbox.checked = true;
      attachPrimerAsFeatureCheckbox.checked = true;
      setPrimerStatus("");
      setPrimerAnalysisPlaceholders();
      savePrimerFromWindowBtn.disabled = false;
      primerSelectionMenu.style.left = `${Math.max(16, Math.min(window.innerWidth - 432, x))}px`;
      primerSelectionMenu.style.top = `${Math.max(16, Math.min(window.innerHeight - 420, y))}px`;
      primerSelectionMenu.style.transform = "none";
      primerSelectionMenu.classList.remove("hidden");
      void analyzePrimerCandidate(candidate);
    }

    async function deleteSelectedPrimerFeature(deletePrimer) {
      const record = getRecord();
      const feature = record && state.selectedFeatureIndex !== null
        ? record.features[state.selectedFeatureIndex]
        : null;
      if (!feature || feature.source !== "user" || !Number(feature.feature_id)) {
        return;
      }

      state.primerDeleteSubmitting = true;
      render();
      mapSelectionSummary.textContent = deletePrimer
        ? "Deleting primer from sequence file and oligostore..."
        : "Removing primer from sequence file...";

      try {
        const response = await fetch(deletePrimerUrl, {
          method: "POST",
          headers: {
            "Content-Type": "application/json",
            "X-CSRFToken": getCsrfToken(),
          },
          body: JSON.stringify({
            feature_id: feature.feature_id,
            delete_primer: deletePrimer,
          }),
        });
        const payload = await response.json().catch(() => ({}));
        if (!response.ok) {
          throw new Error(payload.error || `HTTP ${response.status}`);
        }

        recordDetails.delete(state.recordIndex);
        state.selectedFeatureIndex = null;
        state.pendingFocusFeatureIndex = null;
        state.mapSelectedFeatureIndexes = new Set();
        state.selectedForwardPrimerIndex = null;
        state.selectedReversePrimerIndex = null;
        state.showMapSelectedOnly = false;
        state.tablePage = 1;
        render();

        const needed = getActiveRegionBounds();
        void ensureRecordRegionLoaded(state.recordIndex, needed.regionStart, needed.regionEnd);
      } catch (error) {
        mapSelectionSummary.textContent = deletePrimer
          ? `Could not delete primer: ${error.message}`
          : `Could not remove primer: ${error.message}`;
      } finally {
        state.primerDeleteSubmitting = false;
        render();
      }
    }

    function renderSequenceWindow(record, start, end) {
      const visible = getSequenceSlice(record, start, end);
      if (visible === null) {
        restrictionSiteSummary.textContent = state.loadingRecord
          ? "Loading sequence window..."
          : "Sequence window is outside the loaded region. Adjust pan/zoom to load it.";
        sequenceWindow.textContent = state.loadingRecord ? "Loading sequence..." : "Region not loaded.";
        return;
      }
      const lines = [];
      const hasRestrictionSelection = state.selectedRestrictionEnzymes.size > 0;
      const windowRestrictionHits = getRestrictionHitsInWindow(record, start, end)
        .filter((site) => state.selectedRestrictionEnzymes.has(site.enzyme));
      const highlightedPositions = hasRestrictionSelection
        ? getHighlightedRestrictionPositions(windowRestrictionHits, start, end)
        : new Set();
      const lineNumberLayout = getLineNumberLayout(end);

      for (let i = 0; i < visible.length; i += SEQUENCE_WINDOW_LINE_SIZE) {
        const lineStart = start + i;
        const lineEnd = Math.min(end, lineStart + SEQUENCE_WINDOW_LINE_SIZE - 1);
        const lineSegment = visible.slice(i, i + SEQUENCE_WINDOW_LINE_SIZE);
        const forwardLine = renderStrandChunks(lineSegment, lineStart, highlightedPositions, "forward");
        const reverseLine = renderStrandChunks(lineSegment, lineStart, highlightedPositions, "reverse", complementBase);
        const lineNumber = lineStart.toLocaleString().padStart(lineNumberLayout.width, " ");

        lines.push(`${lineNumber}  5' ${forwardLine} 3'`);
        lines.push(`${lineNumberLayout.prefix}<span class="bottom-strand-row">3' ${reverseLine} 5'</span>`);
        if (hasRestrictionSelection) {
          const annotationRows = buildRestrictionAnnotationRows(windowRestrictionHits, lineStart, lineEnd);
          if (annotationRows) {
            lines.push(
              `${lineNumberLayout.prefix}   <span class="restriction-marker-row"><span class="restriction-cut-marker">${annotationRows.markerText}</span></span>`
            );
            lines.push(
              `${lineNumberLayout.prefix}   <span class="restriction-marker-row">${escapeHtml(annotationRows.labelText)}</span>`
            );
          }
          const lineSites = renderLineRestrictionSites(windowRestrictionHits, lineStart, lineEnd);
          if (lineSites.length) {
            lines.push(`${lineNumberLayout.prefix}<span class="restriction-row">RE: ${lineSites.join(" | ")}</span>`);
          }
        }
      }
      if (hasRestrictionSelection) {
        restrictionSiteSummary.textContent = `${windowRestrictionHits.length.toLocaleString()} restriction site(s) in current window`;
      } else {
        restrictionSiteSummary.textContent = "Select one or more enzymes from the table to visualize binding regions";
      }
      sequenceWindow.innerHTML = lines.join("\n");
    }

    function getSequenceSlice(record, start, end) {
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
    }

    function renderFeatureTrack(record) {
      featureTrack.innerHTML = "";
      forwardPrimerTrack.innerHTML = "";
      reversePrimerTrack.innerHTML = "";
      const colors = [
        "rgba(14, 165, 233, 0.55)",
        "rgba(234, 179, 8, 0.55)",
        "rgba(16, 185, 129, 0.55)",
        "rgba(244, 114, 182, 0.55)",
        "rgba(59, 130, 246, 0.55)",
      ];
      const visibleFeatureIndexes = new Set(getFilteredFeatureIndexes(record));
      const mapStart = state.mapStart;
      const mapEnd = Math.min(record.length, mapStart + state.mapWindowSize - 1);
      const mapLength = mapEnd - mapStart + 1;
      const nonPrimerLanes = [];
      let forwardPrimerCount = 0;
      let reversePrimerCount = 0;

      function buildPrimerArrow(color, isForward) {
        const marker = document.createElement("div");
        marker.style.position = "absolute";
        marker.style.top = "0";
        marker.style.height = "100%";
        marker.style.borderRadius = "9999px";
        marker.style.backgroundColor = color;
        marker.style.border = "1px solid rgba(15, 23, 42, 0.35)";
        marker.style.clipPath = isForward
          ? "polygon(0 0, calc(100% - 8px) 0, 100% 50%, calc(100% - 8px) 100%, 0 100%)"
          : "polygon(8px 0, 100% 0, 100% 100%, 8px 100%, 0 50%)";
        return marker;
      }

      function handleFeatureMouseDown(event, index, feature) {
        if (event.button !== 0) {
          return;
        }
        event.preventDefault();
        applyFeatureSelection(record, index, feature, event.shiftKey);
      }

      const visibleFeatures = record.features
        .map((feature, index) => ({ feature, index }))
        .filter(({ feature }) => {
          const bounds = featureBounds(feature);
          return !(bounds.end < mapStart || bounds.start > mapEnd);
        })
        .sort((a, b) => {
          const aBounds = featureBounds(a.feature);
          const bBounds = featureBounds(b.feature);
          const startDelta = aBounds.start - bBounds.start;
          if (startDelta !== 0) {
            return startDelta;
          }
          return (aBounds.end - aBounds.start) - (bBounds.end - bBounds.start);
        });

      visibleFeatures.forEach(({ feature, index }) => {
        const bounds = featureBounds(feature);
        if (bounds.end < mapStart || bounds.start > mapEnd) {
          return;
        }
        const visibleStart = Math.max(bounds.start, mapStart);
        const visibleEnd = Math.min(bounds.end, mapEnd);
        const left = ((visibleStart - mapStart) / mapLength) * 100;
        const width = (visibleEnd - visibleStart + 1) / mapLength * 100;
        const isSelected = state.selectedFeatureIndex === index;
        const isVisible = visibleFeatureIndexes.has(index);
        if (isPrimerBindingFeature(feature)) {
          const isForward = normalizedStrandValue(feature.strand) !== -1;
          const primerColor = isForward ? "rgba(37, 99, 235, 0.92)" : "rgba(220, 38, 38, 0.92)";
          const primerMarker = buildPrimerArrow(primerColor, isForward);
          primerMarker.className = `cursor-pointer ${isSelected ? "ring-2 ring-primary ring-offset-1 ring-offset-base-200 rounded" : ""} ${isVisible ? "opacity-100" : "opacity-45"}`;
          primerMarker.style.left = `${left}%`;
          primerMarker.style.width = `${Math.max(0.1, Math.min(width, 100 - left))}%`;
          primerMarker.style.minWidth = "12px";
          primerMarker.title = `${feature.label} (${feature.type}) ${feature.start}-${feature.end}`;
          primerMarker.addEventListener("mousedown", (event) => handleFeatureMouseDown(event, index, feature));
          if (isForward) {
            forwardPrimerTrack.appendChild(primerMarker);
            forwardPrimerCount += 1;
          } else {
            reversePrimerTrack.appendChild(primerMarker);
            reversePrimerCount += 1;
          }
          return;
        }

        let lane = 0;
        while (lane < nonPrimerLanes.length && visibleStart <= nonPrimerLanes[lane]) {
          lane += 1;
        }
        nonPrimerLanes[lane] = visibleEnd;

        const marker = document.createElement("div");
        marker.className = `absolute h-6 rounded cursor-pointer ${isSelected ? "ring-2 ring-primary ring-offset-1 ring-offset-base-200" : ""} ${isVisible ? "opacity-95" : "opacity-20"}`;
        marker.style.left = `${left}%`;
        marker.style.width = `${Math.min(width, 100 - left)}%`;
        marker.style.top = `${lane * 22 + 2}px`;
        marker.style.backgroundColor = colors[index % colors.length];
        marker.style.border = "1px solid rgba(15, 23, 42, 0.22)";
        marker.title = `${feature.label} (${feature.type}) ${feature.start}-${feature.end}`;
        marker.addEventListener("mousedown", (event) => handleFeatureMouseDown(event, index, feature));
        featureTrack.appendChild(marker);
      });

      const laneCount = Math.max(1, nonPrimerLanes.length);
      const laneHeight = 22;
      const containerHeight = Math.max(44, laneCount * laneHeight + 6);
      featureTrackContainer.style.height = `${containerHeight}px`;
      primerCountSummary.textContent = `Primers in map range: forward ${forwardPrimerCount.toLocaleString()} | reverse ${reversePrimerCount.toLocaleString()}`;
    }

    function renderWindowFeatureOverlay(record, start, end) {
      windowFeatureOverlay.innerHTML = "";
      windowFeatureOverlay.classList.add("window-feature-overlay-grid");

      if (!state.showWindowFeatureOverlay) {
        windowFeatureOverlaySection.classList.add("hidden");
        toggleWindowFeatureOverlayBtn.textContent = "Show feature lines";
        setWindowFeatureHoverText(null);
        return;
      }

      windowFeatureOverlaySection.classList.remove("hidden");
      toggleWindowFeatureOverlayBtn.textContent = "Hide feature lines";

      const windowLength = Math.max(1, end - start + 1);
      const visibleFeatureIndexes = new Set(getFilteredFeatureIndexes(record));
      const overlappingFeatures = record.features
        .map((feature, index) => ({ feature, index }))
        .filter(({ feature }) => {
          const bounds = featureBounds(feature);
          return !(bounds.end < start || bounds.start > end);
        })
        .sort((a, b) => {
          const aBounds = featureBounds(a.feature);
          const bBounds = featureBounds(b.feature);
          const startDelta = aBounds.start - bBounds.start;
          if (startDelta !== 0) {
            return startDelta;
          }
          return (aBounds.end - aBounds.start) - (bBounds.end - bBounds.start);
        });

      windowFeatureOverlayCount.textContent = `${overlappingFeatures.length.toLocaleString()} in view`;
      if (!overlappingFeatures.length) {
        windowFeatureOverlay.style.height = "2rem";
        const emptyState = document.createElement("div");
        emptyState.className = "absolute inset-0 flex items-center justify-center text-[11px] opacity-60";
        emptyState.textContent = "No annotated features overlap the current sequence window.";
        windowFeatureOverlay.appendChild(emptyState);
        setWindowFeatureHoverText(null);
        return;
      }

      const lanes = [];
      const laneHeight = 10;
      const lineTopOffset = 8;

      overlappingFeatures.forEach(({ feature, index }) => {
        const bounds = featureBounds(feature);
        const visibleStart = Math.max(bounds.start, start);
        const visibleEnd = Math.min(bounds.end, end);
        let lane = 0;
        while (lane < lanes.length && visibleStart <= lanes[lane]) {
          lane += 1;
        }
        lanes[lane] = visibleEnd;

        const left = ((visibleStart - start) / windowLength) * 100;
        const width = ((visibleEnd - visibleStart + 1) / windowLength) * 100;
        const selected = state.selectedFeatureIndex === index;
        const visibleInFilter = visibleFeatureIndexes.has(index);
        const color = isPrimerBindingFeature(feature)
          ? (normalizedStrandValue(feature.strand) === -1 ? "rgba(220, 38, 38, 0.95)" : "rgba(37, 99, 235, 0.95)")
          : "rgba(15, 118, 110, 0.92)";

        const marker = document.createElement("button");
        marker.type = "button";
        marker.className = "absolute rounded-full";
        marker.style.left = `${left}%`;
        marker.style.width = `${Math.max(0.6, Math.min(width, 100 - left))}%`;
        marker.style.minWidth = "6px";
        marker.style.top = `${lane * laneHeight + lineTopOffset}px`;
        marker.style.height = "4px";
        marker.style.backgroundColor = color;
        marker.style.border = "none";
        marker.style.padding = "0";
        marker.style.opacity = visibleInFilter ? "1" : "0.35";
        marker.style.boxShadow = selected
          ? "0 0 0 2px rgba(99, 102, 241, 0.35)"
          : "0 0 0 1px rgba(15, 23, 42, 0.18)";
        marker.style.cursor = "pointer";
        marker.title = `${feature.label} | ${featureDescription(feature)}`;
        marker.addEventListener("mouseenter", () => {
          setWindowFeatureHoverText(feature);
        });
        marker.addEventListener("focus", () => {
          setWindowFeatureHoverText(feature);
        });
        marker.addEventListener("mouseleave", () => {
          const selectedFeature = state.selectedFeatureIndex !== null
            ? record.features[state.selectedFeatureIndex]
            : null;
          setWindowFeatureHoverText(selectedFeature || null);
        });
        marker.addEventListener("blur", () => {
          const selectedFeature = state.selectedFeatureIndex !== null
            ? record.features[state.selectedFeatureIndex]
            : null;
          setWindowFeatureHoverText(selectedFeature || null);
        });
        marker.addEventListener("click", (event) => {
          applyFeatureSelection(record, index, feature, event.shiftKey);
        });
        windowFeatureOverlay.appendChild(marker);
      });

      windowFeatureOverlay.style.height = `${Math.max(32, lanes.length * laneHeight + 16)}px`;
      const selectedFeature = state.selectedFeatureIndex !== null
        ? record.features[state.selectedFeatureIndex]
        : null;
      setWindowFeatureHoverText(selectedFeature || null);
    }

    function renderMapAxis(record) {
      mapTickTrack.innerHTML = "";
      const mapStart = state.mapStart;
      const mapEnd = Math.min(record.length, mapStart + state.mapWindowSize - 1);
      const mapLength = mapEnd - mapStart + 1;
      const targetTicks = 6;
      const step = Math.max(1, Math.floor(mapLength / targetTicks));

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
        mapTickTrack.appendChild(tick);

        const label = document.createElement("div");
        label.className = "absolute top-0 text-[10px] opacity-60";
        label.style.left = `${Math.min(100, Math.max(0, pct))}%`;
        label.style.transform = "translateX(-50%)";
        label.textContent = value.toLocaleString();
        mapTickTrack.appendChild(label);
      });

      mapLeftLabel.textContent = `${mapStart.toLocaleString()} bp`;
      mapRightLabel.textContent = `${mapEnd.toLocaleString()} bp`;
    }

    function renderMapSelectionSummary(record) {
      if (state.selectedFeatureIndex === null || !record.features[state.selectedFeatureIndex]) {
        mapSelectionSummary.textContent = "No feature selected on map.";
        mapSelectionActions.classList.add("hidden");
        return;
      }
      const feature = record.features[state.selectedFeatureIndex];
      mapSelectionSummary.textContent = `Selected: ${feature.label} (${feature.type}) ${feature.start.toLocaleString()}-${feature.end.toLocaleString()} | ${featureLength(feature).toLocaleString()} bp | ${strandLabel(feature.strand)}`;
      const canRemovePrimer = isPrimerBindingFeature(feature) && feature.source === "user" && Number(feature.feature_id) > 0;
      if (!canRemovePrimer) {
        mapSelectionActions.classList.add("hidden");
        return;
      }
      mapSelectionActions.classList.remove("hidden");
      removePrimerFeatureBtn.disabled = state.primerDeleteSubmitting;
      deletePrimerEverywhereBtn.disabled = state.primerDeleteSubmitting || !Number(feature.primer_id);
    }

    function renderFeatureViewer(record) {
      const tableIndexes = getTableFeatureIndexes(record);
      const selectedCount = state.mapSelectedFeatureIndexes.size;
      const modeLabel = state.showMapSelectedOnly ? "map-selected mode" : "all features mode";
      featureCount.textContent = `${tableIndexes.length.toLocaleString()} of ${record.features.length.toLocaleString()} features shown | map-selected: ${selectedCount.toLocaleString()} | ${modeLabel}`;
      if (!tableIndexes.length) {
        const emptyMessage = state.showMapSelectedOnly
          ? "No map-selected features match the current filter."
          : "No features match the current filter.";
        featureTableBody.innerHTML = `<tr><td colspan="6" class="text-center opacity-70 py-3">${emptyMessage}</td></tr>`;
        featurePageInfo.textContent = "Page 0 of 0";
        featurePrevPageBtn.disabled = true;
        featureNextPageBtn.disabled = true;
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
      const pageStart = (state.tablePage - 1) * state.tablePageSize;
      const pageEnd = pageStart + state.tablePageSize;
      const pageIndexes = tableIndexes.slice(pageStart, pageEnd);

      featurePageInfo.textContent = `Page ${state.tablePage.toLocaleString()} of ${totalPages.toLocaleString()}`;
      featurePrevPageBtn.disabled = state.tablePage <= 1;
      featureNextPageBtn.disabled = state.tablePage >= totalPages;

      featureTableBody.innerHTML = "";
      pageIndexes.forEach((index) => {
        const feature = record.features[index];
        const row = document.createElement("tr");
        row.className = `cursor-pointer ${state.selectedFeatureIndex === index ? "bg-base-200" : ""}`;
        row.innerHTML = `
          <td class="max-w-56 truncate" title="${feature.label}">${feature.label}</td>
          <td>${feature.type}</td>
          <td class="font-mono">${feature.start.toLocaleString()}</td>
          <td class="font-mono">${feature.end.toLocaleString()}</td>
          <td class="font-mono">${featureLength(feature).toLocaleString()} bp</td>
          <td>${strandLabel(feature.strand)}</td>
        `;
        row.addEventListener("click", () => {
          state.selectedFeatureIndex = index;
          state.pendingFocusFeatureIndex = index;
          jumpToFeature(record, feature);
          render();
        });
        featureTableBody.appendChild(row);
      });
    }

    function render() {
      const summary = getRecordSummary();
      const record = getRecord();
      clampWindowToRecord();
      clampMapToRecord();
      const recordLength = getCurrentRecordLength();

      windowSizeSlider.min = String(Math.min(MIN_WINDOW_BP, recordLength));
      windowSizeSlider.max = String(Math.min(5000, Math.max(MIN_WINDOW_BP, recordLength)));
      windowSizeSlider.value = String(state.windowSize);

      startSlider.max = String(Math.max(1, recordLength - state.windowSize + 1));
      startSlider.value = String(state.start);

      const end = Math.min(recordLength, state.start + state.windowSize - 1);
      const mapEnd = Math.min(recordLength, state.mapStart + state.mapWindowSize - 1);
      const mapLength = mapEnd - state.mapStart + 1;

      mapStartSlider.max = String(Math.max(1, recordLength - state.mapWindowSize + 1));
      mapStartSlider.value = String(state.mapStart);
      mapRangeLabel.textContent = `Map range: ${state.mapStart.toLocaleString()}-${mapEnd.toLocaleString()} (${mapLength.toLocaleString()} bp)`;

      if (end < state.mapStart || state.start > mapEnd) {
        viewportBox.style.left = "0%";
        viewportBox.style.width = "0%";
      } else {
        const visibleWindowStart = Math.max(state.start, state.mapStart);
        const visibleWindowEnd = Math.min(end, mapEnd);
        const leftPct = ((visibleWindowStart - state.mapStart) / mapLength) * 100;
        const widthPct = ((visibleWindowEnd - visibleWindowStart + 1) / mapLength) * 100;
        viewportBox.style.left = `${leftPct}%`;
        viewportBox.style.width = `${Math.max(widthPct, 0.5)}%`;
      }

      const recordLabel = summary ? summary.id : `Record ${state.recordIndex + 1}`;
      windowLabel.textContent = `${recordLabel} | ${state.start.toLocaleString()}-${end.toLocaleString()} / ${recordLength.toLocaleString()} bp`;
      windowSizeHint.textContent = `${state.windowSize.toLocaleString()} bp visible`;

      if (!record) {
        featureLegend.textContent = state.loadingRecord
          ? "Loading record details..."
          : (state.loadError || "Record details are not available.");
        restrictionSelectionCount.textContent = "0 selected";
        restrictionEnzymeTableBody.innerHTML = `<tr><td colspan="4" class="text-center opacity-70 py-3">${state.loadingRecord ? "Loading..." : "No data loaded."}</td></tr>`;
        restrictionEnzymePageInfo.textContent = "Page 0 of 0";
        restrictionEnzymePrevPageBtn.disabled = true;
        restrictionEnzymeNextPageBtn.disabled = true;
        featureTrack.innerHTML = "";
        forwardPrimerTrack.innerHTML = "";
        reversePrimerTrack.innerHTML = "";
        mapTickTrack.innerHTML = "";
        mapSelectionSummary.textContent = "No feature selected on map.";
        mapSelectionActions.classList.add("hidden");
        primerCountSummary.textContent = "";
        pcrProductSummary.textContent = "";
        pcrProductSequence.textContent = "";
        featureCount.textContent = "";
        featureTableBody.innerHTML = `<tr><td colspan="6" class="text-center opacity-70 py-3">${state.loadingRecord ? "Loading..." : "No feature data loaded."}</td></tr>`;
        featurePageInfo.textContent = "Page 0 of 0";
        featurePrevPageBtn.disabled = true;
        featureNextPageBtn.disabled = true;
        restrictionSiteSummary.textContent = state.loadingRecord ? "Loading sequence window..." : (state.loadError || "Sequence data not loaded.");
        sequenceWindow.textContent = state.loadingRecord ? "Loading sequence..." : "No sequence data loaded.";
        windowFeatureOverlay.innerHTML = "";
        windowFeatureOverlayCount.textContent = "";
        setWindowFeatureHoverText(null);
        savePrimerToOligostoreCheckbox.disabled = true;
        attachPrimerAsFeatureCheckbox.disabled = true;
        reverseComplementPrimerBtn.disabled = true;
        savePrimerFromWindowBtn.disabled = true;
        closePrimerSelectionMenu();
        return;
      }

      const needed = getActiveRegionBounds();
      if (!isRangeCovered(record, needed.regionStart, needed.regionEnd)) {
        void ensureRecordRegionLoaded(state.recordIndex, needed.regionStart, needed.regionEnd);
      }

      featureLegend.textContent = formatFeatures(record.features);
      renderRestrictionEnzymeTable(record);
      renderFeatureTrack(record);
      renderMapAxis(record);
      renderMapSelectionSummary(record);
      renderPCRProduct(record);
      renderFeatureViewer(record);
      renderWindowFeatureOverlay(record, state.start, end);
      renderSequenceWindow(record, state.start, end);
      savePrimerToOligostoreCheckbox.disabled = false;
      attachPrimerAsFeatureCheckbox.disabled = false;
      reverseComplementPrimerBtn.disabled = false;
      savePrimerFromWindowBtn.disabled = state.primerSubmitting || !state.selectedPrimerCandidate;
    }

    recordSummaries.forEach((record, index) => {
      const option = document.createElement("option");
      option.value = String(index);
      option.textContent = `${record.id} (${record.length.toLocaleString()} bp)`;
      recordSelect.appendChild(option);
    });

    recordSelect.addEventListener("change", () => {
      state.recordIndex = Number(recordSelect.value);
      state.windowSize = Math.min(200, getCurrentRecordLength());
      state.start = 1;
      state.mapWindowSize = defaultMapWindowSize(getCurrentRecordLength());
      state.mapStart = 1;
      state.selectedFeatureIndex = null;
      state.tablePage = 1;
      state.showMapSelectedOnly = false;
      state.mapSelectedFeatureIndexes = new Set();
      state.selectedForwardPrimerIndex = null;
      state.selectedReversePrimerIndex = null;
      state.restrictionEnzymeQuery = "";
      state.restrictionTablePage = 1;
      state.selectedRestrictionEnzymes = new Set();
      state.loadError = "";
      state.loadingRecord = false;
      restrictionEnzymeSearchInput.value = state.restrictionEnzymeQuery;
      closePrimerSelectionMenu();
      render();
      const needed = getActiveRegionBounds();
      void ensureRecordRegionLoaded(state.recordIndex, needed.regionStart, needed.regionEnd);
    });

    windowSizeSlider.addEventListener("input", () => {
      closePrimerSelectionMenu();
      state.windowSize = Number(windowSizeSlider.value);
      render();
    });

    startSlider.addEventListener("input", () => {
      closePrimerSelectionMenu();
      state.start = Number(startSlider.value);
      render();
    });

    zoomInBtn.addEventListener("click", () => {
      closePrimerSelectionMenu();
      state.windowSize = Math.max(MIN_WINDOW_BP, Math.floor(state.windowSize / 2));
      render();
    });

    zoomOutBtn.addEventListener("click", () => {
      closePrimerSelectionMenu();
      const maxWindow = Math.min(5000, getCurrentRecordLength());
      state.windowSize = Math.min(maxWindow, state.windowSize * 2);
      render();
    });

    resetBtn.addEventListener("click", () => {
      closePrimerSelectionMenu();
      state.windowSize = Math.min(200, getCurrentRecordLength());
      state.start = 1;
      render();
    });

    mapStartSlider.addEventListener("input", () => {
      closePrimerSelectionMenu();
      state.mapStart = Number(mapStartSlider.value);
      render();
    });

    mapZoomInBtn.addEventListener("click", () => {
      closePrimerSelectionMenu();
      state.mapWindowSize = Math.max(50, Math.floor(state.mapWindowSize / 2));
      render();
    });

    mapZoomOutBtn.addEventListener("click", () => {
      closePrimerSelectionMenu();
      state.mapWindowSize = Math.min(getCurrentRecordLength(), state.mapWindowSize * 2);
      render();
    });

    mapZoomResetBtn.addEventListener("click", () => {
      closePrimerSelectionMenu();
      state.mapWindowSize = defaultMapWindowSize(getCurrentRecordLength());
      state.mapStart = 1;
      render();
    });

    featureSearchInput.addEventListener("input", () => {
      state.featureQuery = featureSearchInput.value;
      state.selectedFeatureIndex = null;
      state.tablePage = 1;
      render();
    });

    featurePageSizeSelect.addEventListener("change", () => {
      state.tablePageSize = Number(featurePageSizeSelect.value);
      state.tablePage = 1;
      render();
    });

    featurePrevPageBtn.addEventListener("click", () => {
      state.tablePage = Math.max(1, state.tablePage - 1);
      render();
    });

    featureNextPageBtn.addEventListener("click", () => {
      state.tablePage = state.tablePage + 1;
      render();
    });

    showMapSelectedBtn.addEventListener("click", () => {
      displayMapSelectedFeaturesInTable();
    });

    showAllFeaturesBtn.addEventListener("click", () => {
      state.showMapSelectedOnly = false;
      state.tablePage = 1;
      render();
    });

    clearMapSelectionBtn.addEventListener("click", () => {
      state.mapSelectedFeatureIndexes = new Set();
      state.selectedFeatureIndex = null;
      state.pendingFocusFeatureIndex = null;
      state.showMapSelectedOnly = false;
      state.selectedForwardPrimerIndex = null;
      state.selectedReversePrimerIndex = null;
      state.tablePage = 1;
      render();
    });

    toggleWindowFeatureOverlayBtn.addEventListener("click", () => {
      state.showWindowFeatureOverlay = !state.showWindowFeatureOverlay;
      render();
    });

    selectTopRestrictionEnzymesBtn.addEventListener("click", () => {
      const record = getRecord();
      if (!record) {
        return;
      }
      applyTopRestrictionSelection(record);
      state.restrictionTablePage = 1;
      render();
    });

    clearRestrictionSelectionBtn.addEventListener("click", () => {
      state.selectedRestrictionEnzymes = new Set();
      state.restrictionTablePage = 1;
      render();
    });

    restrictionEnzymeSearchInput.addEventListener("input", () => {
      state.restrictionEnzymeQuery = restrictionEnzymeSearchInput.value;
      state.restrictionTablePage = 1;
      render();
    });

    restrictionEnzymePrevPageBtn.addEventListener("click", () => {
      state.restrictionTablePage = Math.max(1, state.restrictionTablePage - 1);
      render();
    });

    restrictionEnzymeNextPageBtn.addEventListener("click", () => {
      state.restrictionTablePage = state.restrictionTablePage + 1;
      render();
    });

    closePrimerSelectionMenuBtn.addEventListener("click", () => {
      closePrimerSelectionMenu();
    });

    removePrimerFeatureBtn.addEventListener("click", async () => {
      if (state.primerDeleteSubmitting) {
        return;
      }
      const confirmed = window.confirm("Remove this primer annotation from the sequence file?");
      if (!confirmed) {
        return;
      }
      await deleteSelectedPrimerFeature(false);
    });

    deletePrimerEverywhereBtn.addEventListener("click", async () => {
      if (state.primerDeleteSubmitting) {
        return;
      }
      const confirmed = window.confirm("Delete this primer from the sequence file and oligostore primers?");
      if (!confirmed) {
        return;
      }
      await deleteSelectedPrimerFeature(true);
    });

    document.addEventListener("click", (event) => {
      if (primerSelectionMenu.classList.contains("hidden")) {
        return;
      }
      if (primerSelectionMenu.contains(event.target) || sequenceWindow.contains(event.target)) {
        return;
      }
      closePrimerSelectionMenu();
    });

    document.addEventListener("keydown", (event) => {
      if (event.key === "Escape" && !primerSelectionMenu.classList.contains("hidden")) {
        closePrimerSelectionMenu();
      }
    });

    sequenceWindow.addEventListener("contextmenu", (event) => {
      const candidate = getSelectedPrimerCandidate();
      if (!candidate) {
        return;
      }
      event.preventDefault();
      if (candidate.error) {
        closePrimerSelectionMenu();
        window.alert(candidate.error);
        return;
      }
      openPrimerSelectionMenu(candidate, event.clientX, event.clientY);
    });

    sequenceWindow.addEventListener("mousedown", () => {
      if (!primerSelectionMenu.classList.contains("hidden")) {
        closePrimerSelectionMenu();
      }
    });

    reverseComplementPrimerBtn.addEventListener("click", () => {
      const candidate = state.selectedPrimerCandidate;
      if (!candidate) {
        setPrimerStatus("Select a primer region first.", true);
        return;
      }
      const nextSequence = reverseComplementSequence(candidate.sequence);
      state.selectedPrimerCandidate = {
        ...candidate,
        sequence: nextSequence,
        strand: candidate.strand === -1 ? 1 : -1,
      };
      primerSelectionSummary.textContent = `${candidate.start.toLocaleString()}-${candidate.end.toLocaleString()} | ${nextSequence.length.toLocaleString()} bp | ${strandLabel(state.selectedPrimerCandidate.strand)}`;
      primerSelectionSequence.textContent = nextSequence;
      if (primerNameInput.value.trim() === formatPrimerSelectionName(candidate)) {
        primerNameInput.value = formatPrimerSelectionName(state.selectedPrimerCandidate);
      }
      void analyzePrimerCandidate(state.selectedPrimerCandidate);
      setPrimerStatus("Selection reverse-complemented.");
    });

    savePrimerFromWindowBtn.addEventListener("click", async () => {
      if (state.primerSubmitting || !state.selectedPrimerCandidate) {
        return;
      }
      const candidate = state.selectedPrimerCandidate;
      const primerName = primerNameInput.value.trim();
      const overhangSequence = primerOverhangInput.value.trim();
      const saveToPrimers = savePrimerToOligostoreCheckbox.checked;
      const attachFeature = attachPrimerAsFeatureCheckbox.checked;
      const recordId = getCurrentRecordId();

      if (!primerName) {
        setPrimerStatus("Primer name is required.", true);
        return;
      }
      if (!candidate.sequence) {
        setPrimerStatus("Primer sequence is required.", true);
        return;
      }
      if (candidate.sequence.length > 60) {
        setPrimerStatus("Selected primer is longer than 60 bp and cannot be saved.", true);
        return;
      }
      if (!saveToPrimers && !attachFeature) {
        setPrimerStatus("Select at least one destination.", true);
        return;
      }
      if (attachFeature && !recordId) {
        setPrimerStatus("Could not determine current record for feature attachment.", true);
        return;
      }

      state.primerSubmitting = true;
      render();
      setPrimerStatus("Saving selection...");

      try {
        const response = await fetch(createPrimerUrl, {
          method: "POST",
          headers: {
            "Content-Type": "application/json",
            "X-CSRFToken": getCsrfToken(),
          },
          body: JSON.stringify({
            primer_name: primerName,
            sequence: candidate.sequence,
            overhang_sequence: overhangSequence,
            save_to_primers: saveToPrimers,
            attach_feature: attachFeature,
            record_id: recordId,
            feature_start: candidate.start,
            feature_end: candidate.end,
            feature_strand: candidate.strand,
          }),
        });
        const payload = await response.json().catch(() => ({}));
        if (!response.ok) {
          throw new Error(payload.error || `HTTP ${response.status}`);
        }
        if (payload.attached_feature) {
          recordDetails.delete(state.recordIndex);
          const needed = getActiveRegionBounds();
          void ensureRecordRegionLoaded(state.recordIndex, needed.regionStart, needed.regionEnd);
        }
        const destinations = [];
        if (payload.primer) {
          destinations.push(`oligostore primers as ${payload.primer.name}`);
        }
        if (payload.attached_feature) {
          destinations.push(`sequence file at ${payload.attached_feature.start}-${payload.attached_feature.end}`);
        }
        setPrimerStatus(`Saved to ${destinations.join(" and ")}.`);
      } catch (error) {
        setPrimerStatus(`Could not save selection: ${error.message}`, true);
      } finally {
        state.primerSubmitting = false;
        render();
      }
    });

    state.mapWindowSize = defaultMapWindowSize(getCurrentRecordLength());
    restrictionEnzymeSearchInput.value = state.restrictionEnzymeQuery;
    render();
    const needed = getActiveRegionBounds();
    void ensureRecordRegionLoaded(state.recordIndex, needed.regionStart, needed.regionEnd);

});
