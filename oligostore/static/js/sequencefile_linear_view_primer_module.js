(function () {
  const ns = window.SequenceLinearView = window.SequenceLinearView || {};

  ns.registerPrimerModule = function registerPrimerModule(app) {
    const { state, els, urls } = app;

    app.setPrimerStatus = function setPrimerStatus(message, isError = false) {
      els.primerCreateStatus.textContent = message;
      els.primerCreateStatus.className = `text-xs mt-2 ${isError ? "text-error" : "text-success"}`;
    };

    app.setPCRProductStatus = function setPCRProductStatus(message, isError = false) {
      els.pcrProductSaveStatus.textContent = message;
      els.pcrProductSaveStatus.className = `mt-2 text-xs ${isError ? "text-error" : "text-success"}`;
    };

    app.getPCRProductAutoName = function getPCRProductAutoName(product) {
      const forward = String(product.forwardLabel || "FWD").replace(/\s+/g, "_");
      const reverse = String(product.reverseLabel || "REV").replace(/\s+/g, "_");
      return `${forward}_${reverse}_${product.start}-${product.end}`;
    };

    app.closePrimerSelectionMenu = function closePrimerSelectionMenu() {
      els.primerSelectionMenu.classList.add("hidden");
      state.selectedPrimerCandidate = null;
      state.primerAnalysisLoading = false;
      els.savePrimerFromWindowBtn.disabled = true;
      app.setPrimerStatus("");
    };

    app.setPrimerAnalysisPlaceholders = function setPrimerAnalysisPlaceholders(message = "-") {
      els.primerSelectionTm.textContent = message;
      els.primerSelectionGc.textContent = message;
      els.primerSelectionHairpin.textContent = message;
      els.primerSelectionSelfDimer.textContent = message;
    };

    app.formatPrimerSelectionName = function formatPrimerSelectionName(candidate) {
      const record = app.getRecord();
      const safeRecord = String((record && record.id) || `record_${state.recordIndex + 1}`).replace(/[^A-Za-z0-9_-]/g, "_");
      const direction = candidate.strand === -1 ? "R" : "F";
      return `${safeRecord}_${candidate.start}_${candidate.end}_${direction}`;
    };

    app.getSelectedPrimerCandidate = function getSelectedPrimerCandidate() {
      const selection = window.getSelection();
      if (!selection || selection.rangeCount === 0 || selection.isCollapsed) {
        return null;
      }
      const range = selection.getRangeAt(0);
      if (!els.sequenceWindow.contains(range.commonAncestorContainer)) {
        return null;
      }

      const selectedBaseSpans = Array.from(
        els.sequenceWindow.querySelectorAll('span[data-sequence-base="1"]')
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
        sequence: strand === -1 ? displayedSequence.split("").reverse().join("") : displayedSequence,
      };
    };

    app.analyzePrimerCandidate = async function analyzePrimerCandidate(candidate) {
      state.primerAnalysisLoading = true;
      app.setPrimerAnalysisPlaceholders("Loading...");
      try {
        const response = await fetch(urls.analyzePrimerUrl, {
          method: "POST",
          headers: {
            "Content-Type": "application/x-www-form-urlencoded;charset=UTF-8",
            "X-CSRFToken": app.getCsrfToken(),
          },
          body: new URLSearchParams({ sequence: candidate.sequence }),
        });
        const payload = await response.json().catch(() => ({}));
        if (!response.ok) {
          throw new Error(payload.error || `HTTP ${response.status}`);
        }
        els.primerSelectionTm.textContent = `${payload.tm} deg C`;
        els.primerSelectionGc.textContent = `${Math.round(Number(payload.gc_content || 0) * 100)}%`;
        els.primerSelectionHairpin.textContent = payload.hairpin
          ? `Possible (${payload.hairpin_dg} kcal/mol)`
          : `Not detected (${payload.hairpin_dg} kcal/mol)`;
        els.primerSelectionSelfDimer.textContent = payload.self_dimer
          ? `Possible (${payload.self_dimer_dg} kcal/mol)`
          : `Not detected (${payload.self_dimer_dg} kcal/mol)`;
      } catch (error) {
        app.setPrimerAnalysisPlaceholders("Unavailable");
        app.setPrimerStatus(`Could not analyze selection: ${error.message}`, true);
      } finally {
        state.primerAnalysisLoading = false;
      }
    };

    app.openPrimerSelectionMenu = function openPrimerSelectionMenu(candidate, x, y) {
      state.selectedPrimerCandidate = candidate;
      els.primerSelectionSummary.textContent = `${candidate.start.toLocaleString()}-${candidate.end.toLocaleString()} | ${candidate.sequence.length.toLocaleString()} bp | ${app.strandLabel(candidate.strand)}`;
      els.primerSelectionSequence.textContent = candidate.sequence;
      els.primerNameInput.value = app.formatPrimerSelectionName(candidate);
      els.primerOverhangInput.value = "";
      els.savePrimerToOligostoreCheckbox.checked = true;
      els.attachPrimerAsFeatureCheckbox.checked = true;
      app.setPrimerStatus("");
      app.setPrimerAnalysisPlaceholders();
      els.savePrimerFromWindowBtn.disabled = false;
      els.primerSelectionMenu.style.left = `${Math.max(16, Math.min(window.innerWidth - 432, x))}px`;
      els.primerSelectionMenu.style.top = `${Math.max(16, Math.min(window.innerHeight - 420, y))}px`;
      els.primerSelectionMenu.style.transform = "none";
      els.primerSelectionMenu.classList.remove("hidden");
      void app.analyzePrimerCandidate(candidate);
    };

    app.deleteSelectedPrimerFeature = async function deleteSelectedPrimerFeature(deletePrimer) {
      const record = app.getRecord();
      const feature = record && state.selectedFeatureIndex !== null
        ? record.features[state.selectedFeatureIndex]
        : null;
      if (!feature || feature.source !== "user" || !Number(feature.feature_id)) {
        return;
      }

      state.primerDeleteSubmitting = true;
      app.render();
      els.mapSelectionSummary.textContent = deletePrimer
        ? "Deleting primer from sequence file and oligostore..."
        : "Removing primer from sequence file...";

      try {
        const response = await fetch(urls.deletePrimerUrl, {
          method: "POST",
          headers: {
            "Content-Type": "application/json",
            "X-CSRFToken": app.getCsrfToken(),
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

        app.recordDetails.delete(state.recordIndex);
        state.selectedFeatureIndex = null;
        state.pendingFocusFeatureIndex = null;
        state.mapSelectedFeatureIndexes = new Set();
        state.selectedForwardPrimerIndex = null;
        state.selectedReversePrimerIndex = null;
        state.showMapSelectedOnly = false;
        state.tablePage = 1;
        app.render();

        const needed = app.getActiveRegionBounds();
        void app.ensureRecordRegionLoaded(state.recordIndex, needed.regionStart, needed.regionEnd);
      } catch (error) {
        els.mapSelectionSummary.textContent = deletePrimer
          ? `Could not delete primer: ${error.message}`
          : `Could not remove primer: ${error.message}`;
      } finally {
        state.primerDeleteSubmitting = false;
        app.render();
      }
    };

    app.computePCRProduct = function computePCRProduct(record) {
      if (state.selectedForwardPrimerIndex === null || state.selectedReversePrimerIndex === null) {
        return null;
      }
      const forwardPrimer = record.features[state.selectedForwardPrimerIndex];
      const reversePrimer = record.features[state.selectedReversePrimerIndex];
      if (!forwardPrimer || !reversePrimer) {
        return null;
      }
      const fwdBounds = app.featureBounds(forwardPrimer);
      const revBounds = app.featureBounds(reversePrimer);
      const productStart = fwdBounds.start;
      const productEnd = revBounds.end;
      if (productStart > productEnd) {
        return {
          valid: false,
          reason: "Forward primer is downstream of reverse primer in this linear view.",
        };
      }
      const sequence = app.getSequenceSlice(record, productStart, productEnd);
      if (sequence === null) {
        return {
          valid: false,
          reason: "PCR product sequence is outside the loaded region. Pan/zoom to load that region.",
        };
      }
      return {
        valid: true,
        recordId: app.getCurrentRecordId(),
        start: productStart,
        end: productEnd,
        length: sequence.length,
        sequence,
        forwardLabel: forwardPrimer.label,
        reverseLabel: reversePrimer.label,
        forwardPrimerId: Number.isFinite(Number(forwardPrimer.primer_id)) ? Number(forwardPrimer.primer_id) : null,
        reversePrimerId: Number.isFinite(Number(reversePrimer.primer_id)) ? Number(reversePrimer.primer_id) : null,
        forwardFeatureId: Number.isFinite(Number(forwardPrimer.feature_id)) ? Number(forwardPrimer.feature_id) : null,
        reverseFeatureId: Number.isFinite(Number(reversePrimer.feature_id)) ? Number(reversePrimer.feature_id) : null,
      };
    };

    app.renderPCRProduct = function renderPCRProduct(record) {
      const product = app.computePCRProduct(record);
      if (!product) {
        els.pcrProductSummary.textContent = "PCR product: shift-click one forward and one reverse primer to generate a candidate amplicon.";
        els.pcrProductSequence.textContent = "";
        els.savePcrProductBtn.disabled = true;
        state.lastPCRProductAutoName = "";
        return;
      }
      if (!product.valid) {
        els.pcrProductSummary.textContent = `PCR product unavailable: ${product.reason}`;
        els.pcrProductSequence.textContent = "";
        els.savePcrProductBtn.disabled = true;
        state.lastPCRProductAutoName = "";
        return;
      }
      const autoName = app.getPCRProductAutoName(product);
      if (!els.pcrProductNameInput.value || els.pcrProductNameInput.value === state.lastPCRProductAutoName) {
        els.pcrProductNameInput.value = autoName;
      }
      state.lastPCRProductAutoName = autoName;
      els.pcrProductSummary.textContent = `PCR product: ${product.start.toLocaleString()}-${product.end.toLocaleString()} (${product.length.toLocaleString()} bp) | Fwd: ${product.forwardLabel} | Rev: ${product.reverseLabel}`;
      els.pcrProductSequence.textContent = product.sequence;
      els.savePcrProductBtn.disabled = state.pcrProductSubmitting;
    };

    app.focusWindowOnPCRProduct = function focusWindowOnPCRProduct(record) {
      const product = app.computePCRProduct(record);
      if (!product || !product.valid) {
        return false;
      }
      const maxWindow = Math.min(5000, record.length);
      state.windowSize = Math.max(app.constants.MIN_WINDOW_BP, Math.min(maxWindow, product.length));
      state.start = Math.max(1, Math.min(product.start, record.length - state.windowSize + 1));
      const midpoint = Math.floor((product.start + product.end) / 2);
      const centeredMapStart = midpoint - Math.floor(state.mapWindowSize / 2);
      const maxMapStart = Math.max(1, record.length - state.mapWindowSize + 1);
      state.mapStart = Math.max(1, Math.min(centeredMapStart, maxMapStart));
      return true;
    };

    app.bindPrimerEvents = function bindPrimerEvents() {
      els.closePrimerSelectionMenuBtn.addEventListener("click", () => {
        app.closePrimerSelectionMenu();
      });

      els.removePrimerFeatureBtn.addEventListener("click", async () => {
        if (state.primerDeleteSubmitting) {
          return;
        }
        if (window.confirm("Remove this primer annotation from the sequence file?")) {
          await app.deleteSelectedPrimerFeature(false);
        }
      });

      els.deletePrimerEverywhereBtn.addEventListener("click", async () => {
        if (state.primerDeleteSubmitting) {
          return;
        }
        if (window.confirm("Delete this primer from the sequence file and oligostore primers?")) {
          await app.deleteSelectedPrimerFeature(true);
        }
      });

      els.pcrProductNameInput.addEventListener("input", () => {
        app.setPCRProductStatus("");
      });

      els.savePcrProductBtn.addEventListener("click", async () => {
        const record = app.getRecord();
        const product = record ? app.computePCRProduct(record) : null;
        if (!product || !product.valid) {
          app.setPCRProductStatus("Generate a valid PCR product first.", true);
          return;
        }

        state.pcrProductSubmitting = true;
        app.render();
        app.setPCRProductStatus("Saving PCR product...");

        try {
          const response = await fetch(urls.savePcrProductUrl, {
            method: "POST",
            headers: {
              "Content-Type": "application/json",
              "X-CSRFToken": app.getCsrfToken(),
            },
            body: JSON.stringify({
              name: els.pcrProductNameInput.value.trim(),
              record_id: product.recordId,
              start: product.start,
              end: product.end,
              sequence: product.sequence,
              forward_primer_label: product.forwardLabel,
              reverse_primer_label: product.reverseLabel,
              forward_primer_id: product.forwardPrimerId,
              reverse_primer_id: product.reversePrimerId,
              forward_feature_id: product.forwardFeatureId,
              reverse_feature_id: product.reverseFeatureId,
            }),
          });
          const payload = await response.json().catch(() => ({}));
          if (!response.ok) {
            throw new Error(payload.error || `HTTP ${response.status}`);
          }
          if (payload.pcr_product && payload.pcr_product.name) {
            els.pcrProductNameInput.value = payload.pcr_product.name;
            state.lastPCRProductAutoName = payload.pcr_product.name;
          }
          app.setPCRProductStatus(`Saved as ${payload.pcr_product.name}.`);
        } catch (error) {
          app.setPCRProductStatus(`Could not save PCR product: ${error.message}`, true);
        } finally {
          state.pcrProductSubmitting = false;
          app.render();
        }
      });

      document.addEventListener("click", (event) => {
        if (els.primerSelectionMenu.classList.contains("hidden")) {
          return;
        }
        if (els.primerSelectionMenu.contains(event.target) || els.sequenceWindow.contains(event.target)) {
          return;
        }
        app.closePrimerSelectionMenu();
      });

      document.addEventListener("keydown", (event) => {
        if (event.key === "Escape" && !els.primerSelectionMenu.classList.contains("hidden")) {
          app.closePrimerSelectionMenu();
        }
      });

      els.sequenceWindow.addEventListener("contextmenu", (event) => {
        const candidate = app.getSelectedPrimerCandidate();
        if (!candidate) {
          return;
        }
        event.preventDefault();
        if (candidate.error) {
          app.closePrimerSelectionMenu();
          window.alert(candidate.error);
          return;
        }
        app.openPrimerSelectionMenu(candidate, event.clientX, event.clientY);
      });

      els.sequenceWindow.addEventListener("mousedown", () => {
        if (!els.primerSelectionMenu.classList.contains("hidden")) {
          app.closePrimerSelectionMenu();
        }
      });

      els.reverseComplementPrimerBtn.addEventListener("click", () => {
        const candidate = state.selectedPrimerCandidate;
        if (!candidate) {
          app.setPrimerStatus("Select a primer region first.", true);
          return;
        }
        const nextSequence = app.reverseComplementSequence(candidate.sequence);
        state.selectedPrimerCandidate = {
          ...candidate,
          sequence: nextSequence,
          strand: candidate.strand === -1 ? 1 : -1,
        };
        els.primerSelectionSummary.textContent = `${candidate.start.toLocaleString()}-${candidate.end.toLocaleString()} | ${nextSequence.length.toLocaleString()} bp | ${app.strandLabel(state.selectedPrimerCandidate.strand)}`;
        els.primerSelectionSequence.textContent = nextSequence;
        if (els.primerNameInput.value.trim() === app.formatPrimerSelectionName(candidate)) {
          els.primerNameInput.value = app.formatPrimerSelectionName(state.selectedPrimerCandidate);
        }
        void app.analyzePrimerCandidate(state.selectedPrimerCandidate);
        app.setPrimerStatus("Selection reverse-complemented.");
      });

      els.savePrimerFromWindowBtn.addEventListener("click", async () => {
        if (state.primerSubmitting || !state.selectedPrimerCandidate) {
          return;
        }
        const candidate = state.selectedPrimerCandidate;
        const primerName = els.primerNameInput.value.trim();
        const overhangSequence = els.primerOverhangInput.value.trim();
        const saveToPrimers = els.savePrimerToOligostoreCheckbox.checked;
        const attachFeature = els.attachPrimerAsFeatureCheckbox.checked;
        const recordId = app.getCurrentRecordId();

        if (!primerName) {
          app.setPrimerStatus("Primer name is required.", true);
          return;
        }
        if (!candidate.sequence) {
          app.setPrimerStatus("Primer sequence is required.", true);
          return;
        }
        if (candidate.sequence.length > 60) {
          app.setPrimerStatus("Selected primer is longer than 60 bp and cannot be saved.", true);
          return;
        }
        if (!saveToPrimers && !attachFeature) {
          app.setPrimerStatus("Select at least one destination.", true);
          return;
        }
        if (attachFeature && !recordId) {
          app.setPrimerStatus("Could not determine current record for feature attachment.", true);
          return;
        }

        state.primerSubmitting = true;
        app.render();
        app.setPrimerStatus("Saving selection...");

        try {
          const response = await fetch(urls.createPrimerUrl, {
            method: "POST",
            headers: {
              "Content-Type": "application/json",
              "X-CSRFToken": app.getCsrfToken(),
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
            app.recordDetails.delete(state.recordIndex);
            const needed = app.getActiveRegionBounds();
            void app.ensureRecordRegionLoaded(state.recordIndex, needed.regionStart, needed.regionEnd);
          }
          const destinations = [];
          if (payload.primer) {
            destinations.push(`oligostore primers as ${payload.primer.name}`);
          }
          if (payload.attached_feature) {
            destinations.push(`sequence file at ${payload.attached_feature.start}-${payload.attached_feature.end}`);
          }
          app.setPrimerStatus(`Saved to ${destinations.join(" and ")}.`);
        } catch (error) {
          app.setPrimerStatus(`Could not save selection: ${error.message}`, true);
        } finally {
          state.primerSubmitting = false;
          app.render();
        }
      });
    };
  };
})();
