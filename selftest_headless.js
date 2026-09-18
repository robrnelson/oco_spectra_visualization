/*
 * Headless selftest driver for index.html.
 *
 *     gjs selftest_headless.js            # expects data/selftest_reference.json
 *
 * There is no node on this machine, but gjs (SpiderMonkey) is available. This
 * shims just enough of the DOM that the real, unmodified <script> from
 * index.html can boot with ?selftest=1, so the JS forward model is executed as
 * shipped and compared against the NumPy reference written by validate.py.
 *
 * Run `python3 validate.py` first to generate the reference.
 */

const GLib = imports.gi.GLib;
const ByteArray = imports.byteArray;

/* gjs 1.56 (mozjs60) has no globalThis and no Array.prototype.flatMap, and its
   top-level `this` is a module-ish scope rather than the global object -- the
   Function trick reaches the real global so indirect eval can see the shims.
   Browsers need none of this; it exists only to run the shipped code headlessly. */
const G = Function("return this")();
if (!Array.prototype.flatMap) {
  Array.prototype.flatMap = function (fn, thisArg) {
    return this.reduce((acc, v, i, arr) =>
      acc.concat(fn.call(thisArg, v, i, arr)), []);
  };
}

const HERE = (() => {
  const p = imports.system.programInvocationName;
  const dir = GLib.path_get_dirname(p);
  return GLib.path_is_absolute(dir) ? dir
    : GLib.build_filenamev([GLib.get_current_dir(), dir]);
})();

function readFile(path) {
  const [ok, bytes] = GLib.file_get_contents(path);
  if (!ok) throw new Error(`cannot read ${path}`);
  return ByteArray.toString(bytes);
}

/* ---------- minimal DOM / platform shims ---------------------------------- */

const noop = () => {};

function fakeElement(id) {
  const attrs = {};
  const el = {
    id,
    style: {},
    dataset: {},
    classList: { toggle: noop, add: noop, remove: noop, contains: () => false },
    textContent: "", innerHTML: "", className: "", value: "0",
    offsetWidth: 180, offsetHeight: 20,
    clientWidth: 1100, clientHeight: 620,
    addEventListener: noop, removeEventListener: noop,
    appendChild: noop, setPointerCapture: noop,
    setAttribute: (k, v) => { attrs[k] = v; },
    getAttribute: k => (k in attrs ? attrs[k] : null),
    getBoundingClientRect: () => ({ left: 0, top: 0, width: 1100, height: 620 }),
    getContext: () => makeCtx(el),
    toBlob: noop, click: noop,
  };
  el.parentElement = {
    clientWidth: 1100, clientHeight: 620, addEventListener: noop,
  };
  return el;
}

function makeCtx(canvas) {
  // permissive 2D context: every unknown property is a no-op method
  return new Proxy({ canvas }, {
    get(t, p) { return p in t ? t[p] : noop; },
    set(t, p, v) { t[p] = v; return true; },
  });
}

const elements = {};
G.document = {
  documentElement: fakeElement("html"),
  body: fakeElement("body"),
  getElementById: id => (elements[id] || (elements[id] = fakeElement(id))),
  createElement: tag => fakeElement(tag),
  querySelectorAll: () => [],
};
/* gjs already defines a read-only `window`; add what the app reads onto it. */
if (typeof G.window === "undefined") G.window = {};
try {
  G.window.devicePixelRatio = 1;
  if (!G.window.addEventListener) G.window.addEventListener = noop;
} catch (e) { /* read-only: (window.devicePixelRatio||1) still yields 1 */ }
G.getComputedStyle = () => ({ getPropertyValue: () => "#123456" });
G.location = { search: "?selftest=1", origin: "file://", pathname: "/index.html" };
G.history = { replaceState: noop };
G.navigator = {};
G.performance = { now: () => GLib.get_monotonic_time() / 1000 };
G.ResizeObserver = class { observe() {} disconnect() {} };
G.setTimeout = () => 0;
G.requestAnimationFrame = () => 0;

if (typeof G.console === "undefined") {
  G.console = { log: (...a) => print(a.join(" ")), error: (...a) => print(a.join(" ")),
                         warn: (...a) => print(a.join(" ")) };
}

/* URLSearchParams is absent in gjs 1.56 */
G.URLSearchParams = class {
  constructor(q) {
    this._m = new Map();
    for (const kv of String(q || "").replace(/^\?/, "").split("&")) {
      if (!kv) continue;
      const i = kv.indexOf("=");
      const k = decodeURIComponent(i < 0 ? kv : kv.slice(0, i));
      const v = i < 0 ? "" : decodeURIComponent(kv.slice(i + 1));
      this._m.set(k, v);
    }
  }
  get(k) { return this._m.has(k) ? this._m.get(k) : null; }
  has(k) { return this._m.has(k); }
  set(k, v) { this._m.set(k, String(v)); }
  toString() {
    return [...this._m].map(([k, v]) =>
      encodeURIComponent(k) + "=" + encodeURIComponent(v)).join("&");
  }
};

/* fetch backed by the local filesystem */
let fetchCount = 0;
G.fetch = path => {
  const full = GLib.build_filenamev([HERE, path]);
  try {
    const text = readFile(full);
    fetchCount++;
    return Promise.resolve({ ok: true, status: 200, json: () => Promise.resolve(JSON.parse(text)) });
  } catch (e) {
    return Promise.resolve({ ok: false, status: 404, json: () => Promise.reject(e) });
  }
};

/* ---------- capture what the app reports --------------------------------- */

let bannerText = null;
const bannerEl = document.getElementById("banner");
Object.defineProperty(bannerEl, "textContent", {
  get() { return bannerText || ""; },
  set(v) { bannerText = v; },
  configurable: true,
});

/* ---------- extract and run the real script ------------------------------ */

const html = readFile(GLib.build_filenamev([HERE, "index.html"]));
const m = html.match(/<script>([\s\S]*?)<\/script>\s*<\/body>/);
if (!m) { print("FAIL: could not find the <script> block in index.html"); imports.system.exit(1); }
const source = m[1];

print(`index.html: extracted ${source.length} chars of JS`);

try {
  // indirect eval so the script sees the global shims
  (0, eval)(source);
} catch (e) {
  print("FAIL: script threw while evaluating/booting:");
  print("  " + e);
  if (e.stack) print(e.stack.split("\n").slice(0, 8).map(s => "  " + s).join("\n"));
  imports.system.exit(1);
}

/* The boot path is async; gjs drains the microtask queue after this script
   finishes, so report from an exit handler-ish trailing promise. */
Promise.resolve().then(() => {}).then(() => {}).then(() => {
  // give the app's own promise chain time to settle
}).then(() => {
  GLib.idle_add(GLib.PRIORITY_LOW, () => {
    print(`fetched ${fetchCount} data file(s)`);
    const ms = elements["r_ms"] ? elements["r_ms"].textContent : "";
    if (ms) print(`last reported redraw: ${ms}`);
    if (!bannerText) {
      print("FAIL: no selftest result was reported (did the boot fail?)");
      loop.quit(); failed = true; return GLib.SOURCE_REMOVE;
    }
    print("");
    print(bannerText);
    print("");
    failed = !/^SELFTEST PASS/.test(bannerText);
    print(failed ? "RESULT: FAIL" : "RESULT: PASS");
    loop.quit();
    return GLib.SOURCE_REMOVE;
  });
});

let failed = true;
const loop = GLib.MainLoop.new(null, false);
loop.run();
imports.system.exit(failed ? 1 : 0);
