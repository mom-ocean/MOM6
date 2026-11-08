/* Progressive enhancement: all navigation works without JavaScript. */
// Inner entries are folded first so nested procedure containers stay intact.
for (const container of [
  ...document.querySelectorAll(".api-procedure"),
].reverse()) {
  const target = container.querySelector(
    ':scope > dl[id^="flint-"], :scope > p[id^="flint-"]',
  );
  if (!target) continue;
  const signature =
    target.tagName === "DL" ? target.querySelector(":scope > dt") : target;
  if (!signature) continue;
  const details = document.createElement("details");
  details.className = "procedure";
  const summary = document.createElement("summary");
  summary.innerHTML = signature.innerHTML;
  signature.hidden = true;
  const body = document.createElement("div");
  body.className = "procedure-body";
  body.append(...container.childNodes);
  details.append(summary, body);
  container.replaceWith(details);
}

for (const section of document.querySelectorAll(
  ".module-page main > section",
)) {
  const heading = section.querySelector(":scope > h2");
  if (!heading) continue;
  const details = document.createElement("details");
  details.className = "doc-section";
  details.id = section.id;
  details.open = true;
  const summary = document.createElement("summary");
  summary.append(heading);
  const body = document.createElement("div");
  body.className = "section-body";
  body.append(...section.childNodes);
  details.append(summary, body);
  section.replaceWith(details);
}

function revealTarget(hash = window.location.hash) {
  if (!hash) return;
  let id;
  try {
    id = decodeURIComponent(hash.slice(1));
  } catch {
    return;
  }
  const target = document.getElementById(id);
  if (!target) return;
  for (let parent = target; parent; parent = parent.parentElement) {
    if (parent.tagName === "DETAILS") parent.open = true;
  }
  const visibleTarget = target.closest("details.procedure") || target;
  requestAnimationFrame(() => visibleTarget.scrollIntoView({ block: "start" }));
}
window.addEventListener("hashchange", () => revealTarget());
document.addEventListener("click", (event) => {
  const link = event.target.closest("a[href]");
  if (!link) return;
  const url = new URL(link.href, window.location.href);
  if (
    url.origin === window.location.origin &&
    url.pathname === window.location.pathname &&
    url.search === window.location.search &&
    url.hash
  )
    revealTarget(url.hash);
});

for (const container of document.querySelectorAll(".external-calls")) {
  const details = document.createElement("details");
  details.className = "external-calls";
  const summary = document.createElement("summary");
  summary.textContent = "External calls";
  details.append(summary);
  container.querySelector(".rubric")?.remove();
  details.append(...container.childNodes);
  container.replaceWith(details);
}
revealTarget();

let printClosed = [];
window.addEventListener("beforeprint", () => {
  printClosed = [...document.querySelectorAll("main details:not([open])")];
  for (const details of printClosed) details.open = true;
});
window.addEventListener("afterprint", () => {
  for (const details of printClosed) details.open = false;
  printClosed = [];
});
const filter = document.getElementById("module-filter");
const modules = [...document.querySelectorAll(".module-list li")];
if (window.matchMedia("(max-width: 760px)").matches) {
  document.querySelector(".module-browser").open = false;
  document.querySelector(".guide-browser").open = false;
}
if (filter) {
  document.querySelector(".filter-controls").hidden = false;
  filter.addEventListener("input", () => {
    const query = filter.value.trim().toLowerCase();
    let visible = 0;
    for (const item of modules) {
      item.hidden = !item.textContent.toLowerCase().includes(query);
      if (!item.hidden) visible += 1;
    }
    document.getElementById("filter-status").textContent =
      `${visible} modules found`;
  });
}
