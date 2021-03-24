/*
  Autocomplete code modified from https://www.w3schools.com/howto/howto_js_autocomplete.asp
  P Dexheimer, October 2020
*/

function init_autocomplete(inp, lookahead_url, lookahead_parse_fn) {
    var currentFocus;
    inp.addEventListener("input", function(e) {
        closeAllLists();
        val = this.value;
        if (!val) return false;
        currentFocus = 0;
        var listDiv = document.createElement("div");
        listDiv.setAttribute("id", this.id+"-ac-list");
        listDiv.setAttribute("class", "autocomplete-items");
        this.parentNode.appendChild(listDiv);
        fetch(lookahead_url+"/"+encodeURIComponent(val))
        .then(response => response.json())
        .then(result => {
            for (var i=0; i<result.length; i++) {
                itemDiv = document.createElement("div");
                itemDiv.innerHTML = lookahead_parse_fn(result[i], val);
                itemDiv.addEventListener("click", function(e) {
                    inp.value = this.getElementsByTagName("input")[0].value;
                    closeAllLists();
                });
                listDiv.appendChild(itemDiv);
            }
            setActive(listDiv.getElementsByTagName("div"));
        });
    });
    inp.addEventListener("keydown", function(e) {
        listDiv = document.getElementById(this.id+"-ac-list");
        if (!listDiv) return false;
        items = listDiv.getElementsByTagName("div");
        switch (e.code) {
            case "ArrowDown":
                currentFocus++;
                setActive(items);
                break;
            case "ArrowUp":
                currentFocus--;
                setActive(items);
                break;
            case "Enter":
                if (currentFocus > -1) {
                    e.preventDefault();
                    items[currentFocus].click();
                }
                break;
            case "Escape":
                closeAllLists();
                break;
        }
    });
    function setActive(items) {
        if (!items) {
            currentFocus = -1;
            return false;
        }
        clearActive(items);
        if (currentFocus >= items.length) currentFocus = 0;
        if (currentFocus < 0) currentFocus = items.length - 1;
        items[currentFocus].classList.add("autocomplete-active");
    }
    function clearActive(items) {
        for (var i=0; i<items.length; i++)
            items[i].classList.remove("autocomplete-active");
    }
    function closeAllLists(dontClose) {
        // If dontClose is set, all lists except that one will be closed
        var lists = document.getElementsByClassName("autocomplete-items");
        for (var i=0; i<lists.length; i++) {
            if (dontClose != lists[i]) {
                currentFocus = -1;
                lists[i].parentNode.removeChild(lists[i]);
            }
        }
    }
    document.addEventListener("click", e => closeAllLists(e.target));
}

