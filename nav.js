/* Navigation mobile uniquement — n'interagit pas avec le modèle. */
(function () {
  "use strict";
  const topbar = document.querySelector(".topbar");
  const toggle = document.querySelector(".nav-toggle");
  if (!topbar || !toggle) return;

  function setOpen(open) {
    topbar.classList.toggle("is-open", open);
    document.body.classList.toggle("nav-open", open);
    toggle.setAttribute("aria-expanded", open ? "true" : "false");
  }

  toggle.addEventListener("click", () => {
    setOpen(!topbar.classList.contains("is-open"));
  });

  const backdrop = document.querySelector(".nav-backdrop");
  if (backdrop) backdrop.addEventListener("click", () => setOpen(false));

  topbar.querySelectorAll(".site-links a").forEach((a) => {
    a.addEventListener("click", () => setOpen(false));
  });

  document.addEventListener("keydown", (e) => {
    if (e.key === "Escape") setOpen(false);
  });

  window.addEventListener("resize", () => {
    if (window.matchMedia("(min-width: 900px)").matches) setOpen(false);
  });
})();
