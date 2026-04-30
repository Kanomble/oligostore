(function () {
  const ns = window.SequenceRestriction = window.SequenceRestriction || {};

  const ENZYME_COLORS = [
    "#d97706",
    "#dc2626",
    "#7c3aed",
    "#059669",
    "#2563eb",
    "#c026d3",
    "#0891b2",
    "#65a30d",
  ];

  ns.normalize = function normalize(value) {
    return String(value || "").trim().toLowerCase();
  };

  ns.formatNRun = function formatNRun(count) {
    if (count <= 0) {
      return "";
    }
    if (count <= 8) {
      return "N".repeat(count);
    }
    return `N(${count})`;
  };

  ns.formatCutSite = function formatCutSite(site, cutOffset) {
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
      return `${ns.formatNRun(Math.abs(offset))}/${motif}`;
    }
    return `${motif}/${ns.formatNRun(offset - motif.length)}`;
  };

  ns.getEnzymeStats = function getEnzymeStats(record) {
    const stats = {};
    (record.restriction_sites || []).forEach((site) => {
      const enzyme = String(site.enzyme || "");
      if (!enzyme) {
        return;
      }
      if (!stats[enzyme]) {
        stats[enzyme] = {
          enzyme,
          count: 0,
          site: String(site.site || ""),
          cutOffset: Number(site.cut_offset),
        };
      }
      stats[enzyme].count += 1;
      if (!stats[enzyme].site && site.site) {
        stats[enzyme].site = String(site.site);
      }
      if (!Number.isFinite(stats[enzyme].cutOffset) && Number.isFinite(Number(site.cut_offset))) {
        stats[enzyme].cutOffset = Number(site.cut_offset);
      }
    });
    return stats;
  };

  ns.getSortedEnzymes = function getSortedEnzymes(record) {
    return Object.values(ns.getEnzymeStats(record)).sort((left, right) => {
      const countDelta = right.count - left.count;
      if (countDelta !== 0) {
        return countDelta;
      }
      return left.enzyme.localeCompare(right.enzyme);
    });
  };

  ns.applyTopSelection = function applyTopSelection(record, selectedSet, limit = 5) {
    selectedSet.clear();
    ns.getSortedEnzymes(record).slice(0, limit).forEach((item) => selectedSet.add(item.enzyme));
  };

  ns.getEnzymeColor = function getEnzymeColor(enzyme) {
    const text = String(enzyme || "");
    let hash = 0;
    for (let index = 0; index < text.length; index += 1) {
      hash = ((hash << 5) - hash) + text.charCodeAt(index);
      hash |= 0;
    }
    return ENZYME_COLORS[Math.abs(hash) % ENZYME_COLORS.length];
  };

  ns.createColorSwatch = function createColorSwatch(enzyme) {
    const swatch = document.createElement("span");
    swatch.className = "restriction-enzyme-swatch";
    swatch.style.backgroundColor = ns.getEnzymeColor(enzyme);
    swatch.title = `${enzyme} marker color`;
    return swatch;
  };

  ns.siteTitle = function siteTitle(site) {
    const cutOffset = Number(site.cut_offset);
    const cutSite = ns.formatCutSite(site.site, Number.isFinite(cutOffset) ? cutOffset : NaN);
    return `${site.enzyme} ${Number(site.start).toLocaleString()}-${Number(site.end).toLocaleString()} | ${cutSite}`;
  };
})();
