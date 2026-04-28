(function () {
  const ns = window.SequenceLinearView = window.SequenceLinearView || {};

  function getElement(id) {
    return document.getElementById(id);
  }

  function createApp(raw, config) {
    const recordSummaries = JSON.parse(raw.textContent);
    const initialPcrProductRaw = getElement("initial-pcr-product-data");
    const initialPcrProduct = initialPcrProductRaw ? JSON.parse(initialPcrProductRaw.textContent) : null;
    if (!recordSummaries.length) {
      return null;
    }

    let initialRecordIndex = 0;
    if (initialPcrProduct && initialPcrProduct.record_id) {
      const matchedIndex = recordSummaries.findIndex((record) => String(record.id) === String(initialPcrProduct.record_id));
      if (matchedIndex >= 0) {
        initialRecordIndex = matchedIndex;
      }
    }

    const app = {
      recordSummaries,
      recordDetails: new Map(),
      initialPcrProduct,
      initialPcrProductApplied: false,
      urls: {
        recordDataUrl: config.dataset.recordDataUrl,
        createPrimerUrl: config.dataset.createPrimerUrl,
        deletePrimerUrl: config.dataset.deletePrimerUrl,
        savePcrProductUrl: config.dataset.savePcrProductUrl,
        analyzePrimerUrl: config.dataset.analyzePrimerUrl,
      },
      els: {
        recordSelect: getElement("recordSelect"),
        windowSizeSlider: getElement("windowSizeSlider"),
        startSlider: getElement("startSlider"),
        startInput: getElement("startInput"),
        startHint: getElement("startHint"),
        zoomInBtn: getElement("zoomInBtn"),
        zoomOutBtn: getElement("zoomOutBtn"),
        resetBtn: getElement("resetBtn"),
        mapZoomInBtn: getElement("mapZoomInBtn"),
        mapZoomOutBtn: getElement("mapZoomOutBtn"),
        mapZoomResetBtn: getElement("mapZoomResetBtn"),
        mapStartSlider: getElement("mapStartSlider"),
        mapStartInput: getElement("mapStartInput"),
        mapRangeLabel: getElement("mapRangeLabel"),
        mapTickTrack: getElement("mapTickTrack"),
        mapLeftLabel: getElement("mapLeftLabel"),
        mapRightLabel: getElement("mapRightLabel"),
        toggleCdsFeaturesBtn: getElement("toggleCdsFeaturesBtn"),
        togglePrimerFeaturesBtn: getElement("togglePrimerFeaturesBtn"),
        toggleMiscFeaturesBtn: getElement("toggleMiscFeaturesBtn"),
        primerCountSummary: getElement("primerCountSummary"),
        mapSelectionSummary: getElement("mapSelectionSummary"),
        mapSelectionActions: getElement("mapSelectionActions"),
        removePrimerFeatureBtn: getElement("removePrimerFeatureBtn"),
        deletePrimerEverywhereBtn: getElement("deletePrimerEverywhereBtn"),
        pcrProductSummary: getElement("pcrProductSummary"),
        pcrProductSequence: getElement("pcrProductSequence"),
        pcrProductNameInput: getElement("pcrProductNameInput"),
        savePcrProductBtn: getElement("savePcrProductBtn"),
        pcrProductSaveStatus: getElement("pcrProductSaveStatus"),
        windowLabel: getElement("windowLabel"),
        windowSizeHint: getElement("windowSizeHint"),
        viewportBox: getElement("viewportBox"),
        cdsFeatureTrackContainer: getElement("cdsFeatureTrackContainer"),
        featureTrackContainer: getElement("featureTrackContainer"),
        featureTrack: getElement("primerMiscFeatureTrack"),
        cdsFeatureTrack: getElement("cdsFeatureTrack"),
        primerMiscFeatureTrack: getElement("primerMiscFeatureTrack"),
        featureLegend: getElement("featureLegend"),
        featureSearchInput: getElement("featureSearchInput"),
        featureCount: getElement("featureCount"),
        featureTableBody: getElement("featureTableBody"),
        showMapSelectedBtn: getElement("showMapSelectedBtn"),
        showAllFeaturesBtn: getElement("showAllFeaturesBtn"),
        clearMapSelectionBtn: getElement("clearMapSelectionBtn"),
        featurePageSizeSelect: getElement("featurePageSizeSelect"),
        featurePageInfo: getElement("featurePageInfo"),
        featurePrevPageBtn: getElement("featurePrevPageBtn"),
        featureNextPageBtn: getElement("featureNextPageBtn"),
        selectTopRestrictionEnzymesBtn: getElement("selectTopRestrictionEnzymesBtn"),
        clearRestrictionSelectionBtn: getElement("clearRestrictionSelectionBtn"),
        restrictionEnzymeSearchInput: getElement("restrictionEnzymeSearchInput"),
        restrictionEnzymeTableBody: getElement("restrictionEnzymeTableBody"),
        restrictionEnzymePageInfo: getElement("restrictionEnzymePageInfo"),
        restrictionEnzymePrevPageBtn: getElement("restrictionEnzymePrevPageBtn"),
        restrictionEnzymeNextPageBtn: getElement("restrictionEnzymeNextPageBtn"),
        restrictionSelectionCount: getElement("restrictionSelectionCount"),
        restrictionSiteSummary: getElement("restrictionSiteSummary"),
        sequenceWindow: getElement("sequenceWindow"),
        primerSelectionMenu: getElement("primerSelectionMenu"),
        closePrimerSelectionMenuBtn: getElement("closePrimerSelectionMenuBtn"),
        primerSelectionSummary: getElement("primerSelectionSummary"),
        primerSelectionSequence: getElement("primerSelectionSequence"),
        primerSelectionTm: getElement("primerSelectionTm"),
        primerSelectionGc: getElement("primerSelectionGc"),
        primerSelectionHairpin: getElement("primerSelectionHairpin"),
        primerSelectionSelfDimer: getElement("primerSelectionSelfDimer"),
        primerNameInput: getElement("primerNameInput"),
        primerOverhangInput: getElement("primerOverhangInput"),
        savePrimerToOligostoreCheckbox: getElement("savePrimerToOligostoreCheckbox"),
        attachPrimerAsFeatureCheckbox: getElement("attachPrimerAsFeatureCheckbox"),
        reverseComplementPrimerBtn: getElement("reverseComplementPrimerBtn"),
        savePrimerFromWindowBtn: getElement("savePrimerFromWindowBtn"),
        primerCreateStatus: getElement("primerCreateStatus"),
        toggleWindowFeatureOverlayBtn: getElement("toggleWindowFeatureOverlayBtn"),
        windowFeatureOverlaySection: getElement("windowFeatureOverlaySection"),
        windowFeatureOverlay: getElement("windowFeatureOverlay"),
        windowFeatureOverlayCount: getElement("windowFeatureOverlayCount"),
        windowFeatureHoverDescription: getElement("windowFeatureHoverDescription"),
      },
      constants: {
        MIN_WINDOW_BP: 1,
        PAN_SLIDER_STEPS: 1000,
        SEQUENCE_WINDOW_CHUNK_SIZE: 10,
        SEQUENCE_WINDOW_LINE_SIZE: 100,
        COMPLEMENT_BY_BASE: { A: "T", T: "A", C: "G", G: "C" },
      },
      state: {
        recordIndex: initialRecordIndex,
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
        showCdsFeatures: true,
        showPrimerFeatures: true,
        showMiscFeatures: true,
        pcrProductSubmitting: false,
        lastPCRProductAutoName: "",
      },
    };

    app.getRecord = function getRecord(index = app.state.recordIndex) {
      return app.recordDetails.get(index) || null;
    };
    app.getRecordSummary = function getRecordSummary() {
      return app.recordSummaries[app.state.recordIndex] || null;
    };
    app.getCurrentRecordLength = function getCurrentRecordLength() {
      const loaded = app.getRecord();
      if (loaded) {
        return loaded.length;
      }
      const summary = app.getRecordSummary();
      return summary ? summary.length : 1;
    };
    app.getCurrentRecordId = function getCurrentRecordId() {
      const summary = app.getRecordSummary();
      return summary ? String(summary.id) : "";
    };
    app.getActiveRegionBounds = function getActiveRegionBounds() {
      const recordLength = app.getCurrentRecordLength();
      const sequenceEnd = Math.min(recordLength, app.state.start + app.state.windowSize - 1);
      const mapEnd = Math.min(recordLength, app.state.mapStart + app.state.mapWindowSize - 1);
      const regionStart = Math.max(1, Math.min(app.state.start, app.state.mapStart));
      const regionEnd = Math.max(regionStart, Math.max(sequenceEnd, mapEnd));
      return { regionStart, regionEnd };
    };
    app.isRangeCovered = function isRangeCovered(record, neededStart, neededEnd) {
      return Boolean(record && record.region_start <= neededStart && record.region_end >= neededEnd);
    };
    app.ensureRecordRegionLoaded = async function ensureRecordRegionLoaded(index, neededStart, neededEnd) {
      const cached = app.getRecord(index);
      if (app.isRangeCovered(cached, neededStart, neededEnd)) {
        return;
      }
      if (app.state.loadingRecord) {
        app.state.pendingRegionRequest = { index, start: neededStart, end: neededEnd };
        return;
      }
      if (app.state.recordIndex === index) {
        app.state.loadingRecord = true;
        app.state.loadError = "";
        app.render();
      }
      try {
        const response = await fetch(`${app.urls.recordDataUrl}?record_index=${index}&start=${neededStart}&end=${neededEnd}`);
        if (!response.ok) {
          throw new Error(`HTTP ${response.status}`);
        }
        app.recordDetails.set(index, await response.json());
      } catch (error) {
        if (app.state.recordIndex === index) {
          app.state.loadError = `Could not load sequence record data (${error.message}).`;
        }
      } finally {
        if (app.state.recordIndex === index) {
          app.state.loadingRecord = false;
          app.render();
        }
        const pending = app.state.pendingRegionRequest;
        app.state.pendingRegionRequest = null;
        if (pending) {
          void app.ensureRecordRegionLoaded(pending.index, pending.start, pending.end);
        }
      }
    };
    app.clampPosition = function clampPosition(value, maxStart) {
      return Math.max(1, Math.min(Math.round(Number(value) || 1), maxStart));
    };
    app.defaultMapWindowSize = function defaultMapWindowSize(recordLength) {
      return Math.max(50, Math.min(5000, Math.max(1, recordLength)));
    };
    app.clampWindowToRecord = function clampWindowToRecord() {
      const recordLength = app.getCurrentRecordLength();
      const maxWindow = Math.min(5000, Math.max(app.constants.MIN_WINDOW_BP, recordLength));
      app.state.windowSize = Math.min(Math.max(app.constants.MIN_WINDOW_BP, app.state.windowSize), maxWindow);
      const maxStart = Math.max(1, recordLength - app.state.windowSize + 1);
      app.state.start = app.clampPosition(app.state.start, maxStart);
    };
    app.clampMapToRecord = function clampMapToRecord() {
      const recordLength = app.getCurrentRecordLength();
      app.state.mapWindowSize = Math.min(Math.max(50, app.state.mapWindowSize), Math.max(50, recordLength));
      const maxMapStart = Math.max(1, recordLength - app.state.mapWindowSize + 1);
      app.state.mapStart = app.clampPosition(app.state.mapStart, maxMapStart);
    };
    app.getPanSliderValue = function getPanSliderValue(position, maxStart) {
      if (maxStart <= 1) return 1;
      if (maxStart <= app.constants.PAN_SLIDER_STEPS) return app.clampPosition(position, maxStart);
      return Math.round(((app.clampPosition(position, maxStart) - 1) / (maxStart - 1)) * app.constants.PAN_SLIDER_STEPS);
    };
    app.getPositionFromPanSlider = function getPositionFromPanSlider(rawValue, maxStart) {
      if (maxStart <= 1) return 1;
      if (maxStart <= app.constants.PAN_SLIDER_STEPS) return app.clampPosition(rawValue, maxStart);
      const sliderValue = Math.max(0, Math.min(app.constants.PAN_SLIDER_STEPS, Number(rawValue) || 0));
      return 1 + Math.round((sliderValue / app.constants.PAN_SLIDER_STEPS) * (maxStart - 1));
    };
    app.configurePanControl = function configurePanControl(slider, input, position, maxStart) {
      const usesScaledSlider = maxStart > app.constants.PAN_SLIDER_STEPS;
      slider.min = usesScaledSlider ? "0" : "1";
      slider.max = usesScaledSlider ? String(app.constants.PAN_SLIDER_STEPS) : String(maxStart);
      slider.step = "1";
      slider.value = String(app.getPanSliderValue(position, maxStart));
      slider.disabled = maxStart <= 1;
      input.min = "1";
      input.max = String(maxStart);
      input.step = "1";
      input.value = String(app.clampPosition(position, maxStart));
      input.disabled = maxStart <= 1;
    };
    app.formatFeatures = function formatFeatures(features) {
      if (!features.length) return "No annotated features in this record.";
      const cdsCount = features.filter((feature) => app.isCdsFeature && app.isCdsFeature(feature)).length;
      const primerMiscCount = features.filter((feature) => app.isPrimerOrMiscFeature && app.isPrimerOrMiscFeature(feature)).length;
      return `Displayed lanes: CDS ${cdsCount.toLocaleString()} | primer_bind/misc_feature ${primerMiscCount.toLocaleString()}`;
    };
    app.baseSpan = function baseSpan(base, highlightRestriction = false, position = null, strand = "forward") {
      const value = String(base || "").toUpperCase();
      const restrictionClass = highlightRestriction ? " restriction-hit" : "";
      const metadata = Number.isFinite(position) ? ` data-sequence-base="1" data-position="${position}" data-strand="${strand}"` : "";
      const classPrefix = Number.isFinite(position) ? "sequence-base " : "";
      const className = value === "A" ? "base-a" : value === "C" ? "base-c" : value === "G" ? "base-g" : value === "T" ? "base-t" : "base-n";
      return `<span${metadata} class="${classPrefix}${className}${restrictionClass}">${value}</span>`;
    };
    app.complementBase = function complementBase(base) {
      return app.constants.COMPLEMENT_BY_BASE[String(base || "").toUpperCase()] || "N";
    };
    app.escapeHtml = function escapeHtml(value) {
      return String(value).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;").replace(/"/g, "&quot;").replace(/'/g, "&#39;");
    };
    app.getCsrfToken = function getCsrfToken() {
      const match = document.cookie.match(/(?:^|; )csrftoken=([^;]+)/);
      return match ? decodeURIComponent(match[1]) : "";
    };
    app.reverseComplementSequence = function reverseComplementSequence(sequence) {
      return String(sequence || "").toUpperCase().split("").reverse().map((base) => app.complementBase(base)).join("");
    };

    return app;
  }

  function bindCoreEvents(app) {
    const { state, els } = app;
    els.recordSelect.addEventListener("change", () => {
      state.recordIndex = Number(els.recordSelect.value);
      state.windowSize = Math.min(200, app.getCurrentRecordLength());
      state.start = 1;
      state.mapWindowSize = app.defaultMapWindowSize(app.getCurrentRecordLength());
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
      state.lastPCRProductAutoName = "";
      app.setPCRProductStatus("");
      els.pcrProductNameInput.value = "";
      els.restrictionEnzymeSearchInput.value = "";
      app.closePrimerSelectionMenu();
      app.render();
      const needed = app.getActiveRegionBounds();
      void app.ensureRecordRegionLoaded(state.recordIndex, needed.regionStart, needed.regionEnd);
    });

    els.windowSizeSlider.addEventListener("input", () => { app.closePrimerSelectionMenu(); state.windowSize = Number(els.windowSizeSlider.value); app.render(); });
    els.startSlider.addEventListener("input", () => { app.closePrimerSelectionMenu(); state.start = app.getPositionFromPanSlider(els.startSlider.value, Math.max(1, app.getCurrentRecordLength() - state.windowSize + 1)); app.render(); });
    els.startInput.addEventListener("change", () => { app.closePrimerSelectionMenu(); state.start = app.clampPosition(els.startInput.value, Math.max(1, app.getCurrentRecordLength() - state.windowSize + 1)); app.render(); });
    els.zoomInBtn.addEventListener("click", () => { app.closePrimerSelectionMenu(); state.windowSize = Math.max(app.constants.MIN_WINDOW_BP, Math.floor(state.windowSize / 2)); app.render(); });
    els.zoomOutBtn.addEventListener("click", () => { app.closePrimerSelectionMenu(); state.windowSize = Math.min(Math.min(5000, app.getCurrentRecordLength()), state.windowSize * 2); app.render(); });
    els.resetBtn.addEventListener("click", () => { app.closePrimerSelectionMenu(); state.windowSize = Math.min(200, app.getCurrentRecordLength()); state.start = 1; app.render(); });
    els.mapStartSlider.addEventListener("input", () => { app.closePrimerSelectionMenu(); state.mapStart = app.getPositionFromPanSlider(els.mapStartSlider.value, Math.max(1, app.getCurrentRecordLength() - state.mapWindowSize + 1)); app.render(); });
    els.mapStartInput.addEventListener("change", () => { app.closePrimerSelectionMenu(); state.mapStart = app.clampPosition(els.mapStartInput.value, Math.max(1, app.getCurrentRecordLength() - state.mapWindowSize + 1)); app.render(); });
    els.mapZoomInBtn.addEventListener("click", () => { app.closePrimerSelectionMenu(); state.mapWindowSize = Math.max(50, Math.floor(state.mapWindowSize / 2)); app.render(); });
    els.mapZoomOutBtn.addEventListener("click", () => { app.closePrimerSelectionMenu(); state.mapWindowSize = Math.min(app.getCurrentRecordLength(), state.mapWindowSize * 2); app.render(); });
    els.mapZoomResetBtn.addEventListener("click", () => { app.closePrimerSelectionMenu(); state.mapWindowSize = app.defaultMapWindowSize(app.getCurrentRecordLength()); state.mapStart = 1; app.render(); });
    els.featureSearchInput.addEventListener("input", () => { state.featureQuery = els.featureSearchInput.value; state.selectedFeatureIndex = null; state.tablePage = 1; app.render(); });
    els.featurePageSizeSelect.addEventListener("change", () => { state.tablePageSize = Number(els.featurePageSizeSelect.value); state.tablePage = 1; app.render(); });
    els.featurePrevPageBtn.addEventListener("click", () => { state.tablePage = Math.max(1, state.tablePage - 1); app.render(); });
    els.featureNextPageBtn.addEventListener("click", () => { state.tablePage += 1; app.render(); });
    els.showMapSelectedBtn.addEventListener("click", () => app.displayMapSelectedFeaturesInTable());
    els.showAllFeaturesBtn.addEventListener("click", () => { state.showMapSelectedOnly = false; state.tablePage = 1; app.render(); });
    els.clearMapSelectionBtn.addEventListener("click", () => {
      state.mapSelectedFeatureIndexes = new Set();
      state.selectedFeatureIndex = null;
      state.pendingFocusFeatureIndex = null;
      state.showMapSelectedOnly = false;
      state.selectedForwardPrimerIndex = null;
      state.selectedReversePrimerIndex = null;
      state.lastPCRProductAutoName = "";
      app.setPCRProductStatus("");
      els.pcrProductNameInput.value = "";
      state.tablePage = 1;
      app.render();
    });
    els.toggleWindowFeatureOverlayBtn.addEventListener("click", () => { state.showWindowFeatureOverlay = !state.showWindowFeatureOverlay; app.render(); });
    els.toggleCdsFeaturesBtn.addEventListener("click", () => { state.showCdsFeatures = !state.showCdsFeatures; app.render(); });
    els.togglePrimerFeaturesBtn.addEventListener("click", () => { state.showPrimerFeatures = !state.showPrimerFeatures; app.render(); });
    els.toggleMiscFeaturesBtn.addEventListener("click", () => { state.showMiscFeatures = !state.showMiscFeatures; app.render(); });
    els.selectTopRestrictionEnzymesBtn.addEventListener("click", () => { const record = app.getRecord(); if (record) { app.applyTopRestrictionSelection(record); state.restrictionTablePage = 1; app.render(); } });
    els.clearRestrictionSelectionBtn.addEventListener("click", () => { state.selectedRestrictionEnzymes = new Set(); state.restrictionTablePage = 1; app.render(); });
    els.restrictionEnzymeSearchInput.addEventListener("input", () => { state.restrictionEnzymeQuery = els.restrictionEnzymeSearchInput.value; state.restrictionTablePage = 1; app.render(); });
    els.restrictionEnzymePrevPageBtn.addEventListener("click", () => { state.restrictionTablePage = Math.max(1, state.restrictionTablePage - 1); app.render(); });
    els.restrictionEnzymeNextPageBtn.addEventListener("click", () => { state.restrictionTablePage += 1; app.render(); });
  }

  function render(app) {
    const { state, els } = app;
    const summary = app.getRecordSummary();
    const record = app.getRecord();
    app.clampWindowToRecord();
    app.clampMapToRecord();
    const recordLength = app.getCurrentRecordLength();

    els.windowSizeSlider.min = "1";
    els.windowSizeSlider.max = String(Math.min(5000, Math.max(1, recordLength)));
    els.windowSizeSlider.value = String(state.windowSize);
    const maxStart = Math.max(1, recordLength - state.windowSize + 1);
    app.configurePanControl(els.startSlider, els.startInput, state.start, maxStart);

    const end = Math.min(recordLength, state.start + state.windowSize - 1);
    const mapEnd = Math.min(recordLength, state.mapStart + state.mapWindowSize - 1);
    const mapLength = mapEnd - state.mapStart + 1;
    const maxMapStart = Math.max(1, recordLength - state.mapWindowSize + 1);
    app.configurePanControl(els.mapStartSlider, els.mapStartInput, state.mapStart, maxMapStart);
    els.mapRangeLabel.textContent = `Map range: ${state.mapStart.toLocaleString()}-${mapEnd.toLocaleString()} (${mapLength.toLocaleString()} bp)`;

    if (end < state.mapStart || state.start > mapEnd) {
      els.viewportBox.style.left = "0%";
      els.viewportBox.style.width = "0%";
    } else {
      const visibleWindowStart = Math.max(state.start, state.mapStart);
      const visibleWindowEnd = Math.min(end, mapEnd);
      els.viewportBox.style.left = `${((visibleWindowStart - state.mapStart) / mapLength) * 100}%`;
      els.viewportBox.style.width = `${Math.max(((visibleWindowEnd - visibleWindowStart + 1) / mapLength) * 100, 0.5)}%`;
    }

    const recordLabel = summary ? summary.id : `Record ${state.recordIndex + 1}`;
    els.windowLabel.textContent = `${recordLabel} | ${state.start.toLocaleString()}-${end.toLocaleString()} / ${recordLength.toLocaleString()} bp`;
    els.windowSizeHint.textContent = `${state.windowSize.toLocaleString()} bp visible`;
    els.startHint.textContent = maxStart > app.constants.PAN_SLIDER_STEPS
      ? `Start ${state.start.toLocaleString()} of ${maxStart.toLocaleString()} possible positions | slider is scaled for long sequences`
      : `Start ${state.start.toLocaleString()} of ${maxStart.toLocaleString()} possible positions`;

    if (!record) {
      els.featureLegend.textContent = state.loadingRecord ? "Loading record details..." : (state.loadError || "Record details are not available.");
      els.restrictionSelectionCount.textContent = "0 selected";
      els.restrictionEnzymeTableBody.innerHTML = `<tr><td colspan="4" class="text-center opacity-70 py-3">${state.loadingRecord ? "Loading..." : "No data loaded."}</td></tr>`;
      els.restrictionEnzymePageInfo.textContent = "Page 0 of 0";
      els.restrictionEnzymePrevPageBtn.disabled = true;
      els.restrictionEnzymeNextPageBtn.disabled = true;
      els.featureTrack.innerHTML = "";
      els.cdsFeatureTrack.innerHTML = "";
      els.primerMiscFeatureTrack.innerHTML = "";
      els.mapTickTrack.innerHTML = "";
      els.mapSelectionSummary.textContent = "No feature selected on map.";
      els.mapSelectionActions.classList.add("hidden");
      els.primerCountSummary.textContent = "";
      els.pcrProductSummary.textContent = "";
      els.pcrProductSequence.textContent = "";
      els.savePcrProductBtn.disabled = true;
      els.featureCount.textContent = "";
      els.featureTableBody.innerHTML = `<tr><td colspan="6" class="text-center opacity-70 py-3">${state.loadingRecord ? "Loading..." : "No feature data loaded."}</td></tr>`;
      els.featurePageInfo.textContent = "Page 0 of 0";
      els.featurePrevPageBtn.disabled = true;
      els.featureNextPageBtn.disabled = true;
      els.restrictionSiteSummary.textContent = state.loadingRecord ? "Loading sequence window..." : (state.loadError || "Sequence data not loaded.");
      els.sequenceWindow.textContent = state.loadingRecord ? "Loading sequence..." : "No sequence data loaded.";
      els.windowFeatureOverlay.innerHTML = "";
      els.windowFeatureOverlayCount.textContent = "";
      app.setWindowFeatureHoverText(null);
      els.savePrimerToOligostoreCheckbox.disabled = true;
      els.attachPrimerAsFeatureCheckbox.disabled = true;
      els.reverseComplementPrimerBtn.disabled = true;
      els.savePrimerFromWindowBtn.disabled = true;
      app.closePrimerSelectionMenu();
      return;
    }

    const needed = app.getActiveRegionBounds();
    if (!app.isRangeCovered(record, needed.regionStart, needed.regionEnd)) {
      void app.ensureRecordRegionLoaded(state.recordIndex, needed.regionStart, needed.regionEnd);
    }

    if (app.applyInitialPCRProductSelection && !app.initialPcrProductApplied) {
      app.applyInitialPCRProductSelection(record);
    }

    els.featureLegend.textContent = app.formatFeatures(record.features);
    app.renderRestrictionEnzymeTable(record);
    app.renderFeatureTrack(record);
    app.renderMapAxis(record);
    app.renderMapSelectionSummary(record);
    app.renderPCRProduct(record);
    app.renderFeatureViewer(record);
    app.renderWindowFeatureOverlay(record, state.start, end);
    app.renderSequenceWindow(record, state.start, end);
    els.savePrimerToOligostoreCheckbox.disabled = false;
    els.attachPrimerAsFeatureCheckbox.disabled = false;
    els.reverseComplementPrimerBtn.disabled = false;
    els.savePrimerFromWindowBtn.disabled = state.primerSubmitting || !state.selectedPrimerCandidate;
  }

  document.addEventListener("DOMContentLoaded", () => {
    const raw = getElement("sequence-records-data");
    const config = getElement("sequence-linear-config");
    if (!raw || !config) {
      return;
    }

    const app = createApp(raw, config);
    if (!app) {
      return;
    }

    ns.registerFeatureModule(app);
    ns.registerPrimerModule(app);
    app.render = () => render(app);
    bindCoreEvents(app);
    app.bindPrimerEvents();

    app.recordSummaries.forEach((record, index) => {
      const option = document.createElement("option");
      option.value = String(index);
      option.textContent = `${record.id} (${record.length.toLocaleString()} bp)`;
      app.els.recordSelect.appendChild(option);
    });

    app.els.recordSelect.value = String(app.state.recordIndex);
    app.state.mapWindowSize = app.defaultMapWindowSize(app.getCurrentRecordLength());
    app.els.restrictionEnzymeSearchInput.value = app.state.restrictionEnzymeQuery;
    app.render();
    const needed = app.getActiveRegionBounds();
    void app.ensureRecordRegionLoaded(app.state.recordIndex, needed.regionStart, needed.regionEnd);
  });
})();
