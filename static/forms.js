/**
 * Given a string that looks like ``some[text][with][brackets]``, return the
 * part inside the first set of square brackets `[]`
 *
 * @param {String} str - string to check
 * @param {Number} n - index of key to return
 * @returns {String} key
 */
function nth_key(str, n=0) {
    const subscript_re = /\[([^\]]*)\]/g; // match text inside these: []

    // each match looks like ["[text]", "text"], so matches[0][1] gets the
    // first subscript, not including brackets
    const matches = [...str.matchAll(subscript_re)];
    const key = matches[n][1];

    return key;
}

/**
 * Shortcut for :js:func:`nth_key` with ``n=0``.
 */
function first_key(str) {
    return nth_key(str);
}

/**
 * Shortcuts to select jQuery elements by name or other attribute
 */

/**
 * Return jQuery object with name selector; for name-only selections
 *
 * @param {String} sel - name to select
 * @param {String} rel - relational operator
 * @returns {jQuery} elements
 */
function jq_name(sel, rel='=') {
    return $(sel_name(sel, rel));
}

/**
 * Return an name selector string; for use in compound selections
 *
 * @returns {String} CSS selector
 */
function sel_name(sel, rel='=') {
    return `#${active_page} [name${rel}"${sel}"]`;
}

/**
 * Select elements by a single attribute
 *
 * @returns {jQuery} elements
 */
function jq_attr(attr, sel, rel='=') {
    return $(sel_attr(attr, sel, rel));
}

/**
 * Return an attribute selector string
 *
 * @returns {String} CSS selector
 */
function sel_attr(attr, sel, rel='=') {
    return `#${active_page} [${attr}${rel}"${sel}"]`;
}

/**
 * Return an ID selector that is unique across all pages
 *
 * @returns {String} CSS selector
 */
function sel_id(id, rel='=') {
    if (['=', '^='].includes(rel))
        return `[id${rel}"${active_page}_${id}"]`;
    return `[id${rel}"${id}"]`;
}

/**
 * Return a jQuery objected selected by ID (unique across pages)
 *
 * @returns {jQuery} element
 */
function jq_id(id, rel='=') {
    return $(sel_id(id, rel));
}

/**
 * Select an input's label
 *
 * @param {String} for_id - ID of ``<input>`` element
 * @param {String} [rel]
 * @returns {jQuery} element
 */
function jq_label(for_id, rel='=') {
    return $(sel_label(for_id, rel));
}

/**
 * @param {String} for_id - ID of ``<input>`` element
 * @param {String} [rel]
 * @returns {String} CSS selector
 */
function sel_label(for_id, rel='=') {
    return `label[for${rel}"${active_page}_${for_id}"]`;
}

/**
 *  Callback functions should take 2 args:
 *      target  the element to update
 *      src     the form input whose value is used for updating
 */

/**
 * Converts a number to a comma separated list of numbers from 1 to the
 * target's number. E.g., ``4`` becomes ``1,2,3,4``.
 */
function num_to_range_list(target, src) {
    let from = $(target).attr('data-range-start');
    if (!from)
        from = 0;

    const len = Number.parseInt( src.val() ) - from;

    let nums = Array(len);
    for (let i = 0; i < len; ++i)
        nums[i] = from++;

    target.html(nums.join(','));
}

/**
 * Applies updates from shell script template to the code textarea
 */
function update_code() {
    jq_name('code').val( jq_id('script_tpl').text() );
}

/**
 * Returns a callback function to toggle visibility of an email option
 * subsection.
 *
 * @params {String} assign_page - page_id of the currently active page
 * @return {function} callback - the toggler
 */
function toggle_email(assign_page) {
    var assigned_page = assign_page;
    if (!assign_page)
        console.log('Missing page assignment for toggle_email()');

    let box_sel = `#${assigned_page} [name='PBS[mail][use]']`;
    var box = $(box_sel);
    var email_opts = $(`#${assigned_page} .email-options`);

    if (!box.length)
        console.log(`Missing box: ${box_sel}`);

    box = box[0];

    return function() {
        if (box.checked)
            email_opts.show();
        else
            email_opts.hide();
    }
}

var active_page = null;
var loaded_pages = [];

/**
 * Handles page changes, initializes event handlers
 *
 * @param {jQuery.Event} event - this is ignored
 * @param {Object} ui - see jQuery mobile `pagecontainerload`_
 */
function setup_form(event, ui) {
    active_page = ui.toPage[0].id;

    if (loaded_pages.includes(active_page))
        return;

    if (['home', '404'].includes(active_page))
        return;

    $(`#${active_page} [data-visibility]`).each(function(index, elem) {
        let rules = JSON.parse($(elem).attr('data-visibility'));
        console.log(rules);
    });
}

$(document).on('pageload', function(event, ui) {
    active_page = ui.toPage[0].id;
    let index = loaded_pages.indexOf(active_page);
    if (index !== -1)
        loaded_pages.splice(index, 1); // remove it to allow re-loading
});
$(document).on('pagecontainerchange', setup_form);
