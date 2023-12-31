
let DO_PREDICT = false;

function makeElement(what, cssclass, elt) {
  var e = document.createElement(what);
  if (cssclass) e.setAttribute('class', cssclass);
  if (elt) elt.appendChild(e);
  return e;
}
function DIV(cssclass, elt) { return makeElement('DIV', cssclass, elt); }
function CANVAS(cssclass, elt) { return makeElement('CANVAS', cssclass, elt); }
function SPAN(cssclass, elt) { return makeElement('SPAN', cssclass, elt); }
function BR(cssclass, elt) { return makeElement('BR', cssclass, elt); }
function TEXT(contents, elt) {
  var e = document.createTextNode(contents);
  if (elt) elt.appendChild(e);
  return e;
}

function FetchLines(url, k) {
  let xhr = new XMLHttpRequest;
  xhr.open('GET', url);
  xhr.onreadystatechange = e => {
    if (xhr.readyState === XMLHttpRequest.DONE) {
      k(xhr.responseText);
    }
  };
  xhr.send();
}

function IsWord(w) {
  return WORD_RE.test(w);
}

function WordAt(idx) {
  let elt = document.getElementById('word' + idx);
  return elt.value;
}

// PERF Could cache the list...
function SetSimilar(idx, value) {
  let elt = document.getElementById('under' + idx);
  FetchLines('/similar/' + WORD[idx] + '/' + value,
             text => {
               elt.innerHTML = '';
               let lines = text.split(/\r?\n/);
               for (let line of lines) {
                 TEXT(line, elt);
                 BR('', elt);
               }
             });
}

/*
// Perhaps should do this in batch, as the model
// can efficiently predict many at once.
function Predict(idx) {
  // The starting phrase should be at least 7 words,
  // so that we always have enough context. Best if
  // the end is informative.
  let phrase = 'an acronym defining the word ' + WORD + ' is';
  for (let i = 0; i < idx; i++) {
    phrase += ' ' + WordAt(i);
  }
  let elt = document.getElementById('pred' + idx);

  FetchLines('/next/' + phrase,
             text => {
               elt.innerHTML = '';
               let lines = text.split(/\r?\n/);
               for (let line of lines) {
                 TEXT(line, elt);
                 BR('', elt);
               }
             });
}
*/

function Predict(idx) {
  if (!DO_PREDICT) return;

  let phrase = '';
  for (let i = 0; i < WORD.length; i++) {
    if (i != 0) phrase += ' ';
    phrase += WordAt(i);
  }
  let elt = document.getElementById('pred' + idx);

  FetchLines('/guess/' + WORD + '/' + WORD[idx] + '/' + idx + '/' + phrase,
             text => {
               elt.innerHTML = '';
               let lines = text.split(/\r?\n/);
               for (let line of lines) {
                 TEXT(line, elt);
                 BR('', elt);
               }
             });
}

function InputChange(idx, elt, value) {
  // XXX check start char
  if (IsWord(value)) {
    elt.style.background = '#DFD';
    elt.style.border = '2px solid #090';
    SetSimilar(idx, value);
    // All predictions are affected.
    for (let i = 0; i < WORD.length; i++) {
      Predict(i);
    }
  } else {
    elt.style.background = '#FDD';
    elt.style.border = '2px solid #900';
  }
}

function Fill() {
  for (let idx = 0; idx < WORD.length; idx++) {
    let elt = document.getElementById('word' + idx);
    elt.value = WORD;
    SetSimilar(idx, WORD);
  }
  for (let i = 0; i < WORD.length; i++) {
    Predict(i);
  }
}

function Start() {
  // Register handlers
  for (let idx = 0; idx < WORD.length; idx++) {
    let elt = document.getElementById('word' + idx);
    if (!elt) console.error('missing elt?');

    elt.addEventListener('input', () => InputChange(idx, elt, elt.value));
    InputChange(idx, elt, '');
  }

  Fill();
}

